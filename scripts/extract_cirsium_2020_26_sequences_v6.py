"""Extract a *bounded, exact* set of real Park 2020 Cirsium unigene sequences.

This bridges publisher-supplied BLAST annotation XLSX rows (already audited
into data/molecular/cirsium_japonicum_2020_blast_candidate_hits_v6.json) to
the publisher's separate Supplementary Information 2 transcript FASTA. A
sequence entry remains a **BLAST-annotated candidate**, NOT a validated DFR,
ANS or pigmentation locus. Neither paper nor these sequences describe the
focal Ryukyu *brevicaule/irumtiense* six-individual RNA libraries.

Source archive: Europe PMC public PMC7675961 supplementaryFiles.
The archive SHA256 and exact unigenes are fixed before extraction; do not
silently accept a changed archive or unrelated FASTA.
"""
import argparse
from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import re
import sys
import urllib.error
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7675961/supplementaryFiles"
REGISTER = ROOT / "data/molecular/cirsium_japonicum_2020_blast_candidate_hits_v6.json"
TEXT_NAME = "ao0c04001_si_002.txt"
MAX_ZIP = 50_000_000
MAX_EXPANDED = 180_000_000
HEADER_PATTERN = re.compile(r"^(TR\d+\|c\d+_g\d+_i\d+)\s+len=(\d+)(?:\s|$)")


def load_registry(path=REGISTER):
    d=json.loads(Path(path).read_text(encoding="utf-8"))
    if d.get("evidence_contract") != "anto_cirsium_japonicum_2020_source_blast_hits_v6":
        raise ValueError("Unrecognized registration contract")
    records=d.get("records")
    if not isinstance(records,list) or len(records) != 26 or d.get("selected_rows") != 26:
        raise ValueError("A frozen set of 26 publisher-annotation candidates is required")
    ids=[r["unigene_id"] for r in records]
    if len(set(ids)) != 26 or any(not re.fullmatch(r"TR\d+\|c\d+_g\d+_i\d+",x) for x in ids):
        raise ValueError("Unregistered/duplicate unigene IDs")
    if not re.fullmatch(r"[a-f0-9]{64}",d["source_zip_sha256"]):
        raise ValueError("Source checksum missing")
    return d


def extract_text_stream(raw, targets, expected_total_headers=104890):
    """One pass through source FASTA; strict match for all declared transcript IDs."""
    found={}
    current=None
    seq=[]
    declared=None
    seen_headers=0
    def emit():
        nonlocal current,seq,declared
        if current is None:
            return
        bases="".join(seq)
        if len(bases)!=declared:
            raise ValueError(f"Publisher FASTA length mismatch for {current}: {len(bases)} != {declared}")
        if not bases or any(ch not in "ACGTUNRYKMSWBDHV-" for ch in bases):
            raise ValueError("Unexpected nucleotide alphabet in selected transcript")
        if current in found:
            raise ValueError("Duplicate target transcript accession in publisher FASTA")
        found[current] = {
            "id":current,
            "sequence":bases,
            "source_fasta_length":declared,
            "sha256_nucleotide":sha256(bases.encode("ascii")).hexdigest(),
        }
    for rawline in raw:
        line=rawline.decode("utf-8", errors="replace").strip()
        if line.startswith(">"):
            if current is not None:
                emit()
            seen_headers += 1
            m=HEADER_PATTERN.match(line[1:])
            if m is None:
                raise ValueError(f"Unparseable publisher unigene FASTA header: {line[:150]}")
            ident,decl_len=m.group(1),int(m.group(2))
            current=ident if ident in targets else None
            declared=decl_len if current else None
            seq=[]
        elif current and line:
            seq.append(line.upper())
    if current:
        emit()
    if seen_headers != expected_total_headers:
        raise ValueError(f"Source unigene count drift: {seen_headers}")
    if set(found) != set(targets):
        raise ValueError(f"Missing {len(set(targets)-set(found))} target transcripts: {sorted(set(targets)-set(found))}")
    return found,seen_headers


def extract_archive(data, registry):
    if len(data)>MAX_ZIP or not zipfile.is_zipfile(BytesIO(data)):
        raise ValueError("Archive exceeds budget or is not a ZIP")
    actual_sha=sha256(data).hexdigest()
    if actual_sha != registry["source_zip_sha256"]:
        raise ValueError(f"Source ZIP checksum mismatch; requires independent audit: {actual_sha}")
    with zipfile.ZipFile(BytesIO(data),"r") as z:
        members=z.infolist()
        if sum(m.file_size for m in members)>MAX_EXPANDED:
            raise ValueError("Expanded archive exceeds budget")
        target=next((m for m in members if Path(m.filename).name==TEXT_NAME),None)
        if target is None:
            raise ValueError("Frozen publisher supplementary TXT missing")
        with z.open(target) as stream:
            found,total_headers=extract_text_stream(
                stream,{r["unigene_id"] for r in registry["records"]})
    rows=[]
    for r in registry["records"]:
        c=found[r["unigene_id"]]
        rows.append({
            **c,
            "published_blast_annotation":r["annotated_best_blast_hit"],
            "published_blast_subject":r["blast_subject_identifier"],
            "published_blast_evalue":r["blast_evalue_source_string"],
            "published_annotation_sheet_row":r["supplement_sheet_row"],
            "family_by_publisher_annotation_only":r["annotated_family_candidate"],
            "orthology_verified":False,
            "coding_region_resolved":False,
            "protein_function_verified":False,
            "flower_expression_measured":False,
            "Ryukyu_2026_focal_sample":False,
        })
    return {
        "contract":"anto_park2020_26_source_transcripts_v6",
        "source_study":"10.1021/acsomega.0c04001",
        "source_archive":SOURCE_URL,
        "source_zip_sha256":actual_sha,
        "source_fasta_entry_count":total_headers,
        "selected_annotation_rows":len(rows),
        "sequence_extracted_rows":len(rows),
        "selected_total_nucleotide_chars":sum(r["source_fasta_length"] for r in rows),
        "distinct_sequences":len({r["sha256_nucleotide"] for r in rows}),
        "sample_id_origin":"2020 Cirsium japonicum assembled transcriptome, NOT 2026 Ryukyu six runs",
        "state":"PUBLISHER_FASTA_SEQUENCE_EXTRACTED_ORTHOLOGY_UNVERIFIED",
        "focal_Ryukyu_reexpression_confirmed":False,
        "records":rows,
    }


def fetch_zip():
    req=urllib.request.Request(SOURCE_URL,headers={
        "User-Agent":"anto-reproducible-2020-source-check/1.0",
        "Accept":"application/zip,application/octet-stream"})
    with urllib.request.urlopen(req,timeout=50) as response:
        declared=response.headers.get("Content-Length")
        if declared and int(declared)>MAX_ZIP:
            raise ValueError("Server-declared archive size too large")
        data=response.read(MAX_ZIP+1)
        if len(data)>MAX_ZIP:
            raise ValueError("Archive download exceeds cap")
        return data


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fetch-public",action="store_true")
    p.add_argument("--archive",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--print-transfer-json",action="store_true",
                   help="Print one source-bounded JSON payload between sentinels for audited GitHub Connector ingestion")
    args=p.parse_args()
    if args.fetch_public == (args.archive is not None):
        raise ValueError("Choose exactly one input mode")
    reg=load_registry()
    data=fetch_zip() if args.fetch_public else args.archive.read_bytes()
    outcome=extract_archive(data,reg)
    compact=json.dumps(outcome,separators=(",",":"),ensure_ascii=False)
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(outcome,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("Extracted 26 published transcript FASTA records; biological orthology remains unverified.")
    print(f"Source SHA256={outcome['source_zip_sha256']}; distinct sequences={outcome['distinct_sequences']}; total_nt={outcome['selected_total_nucleotide_chars']}")
    if args.print_transfer_json:
        print("ANTO_CIRSIUM_2020_TRANSFER_BEGIN")
        without_records = {k:v for k,v in outcome.items() if k != "records"}
        print("ANTO_SOURCE_SUMMARY_JSON:"+json.dumps(without_records,separators=(",",":"),ensure_ascii=False))
        for record in outcome["records"]:
            print("ANTO_SOURCE_RECORD_JSON:"+json.dumps(record,separators=(",",":"),ensure_ascii=False))
        print("ANTO_CIRSIUM_2020_TRANSFER_END")


if __name__=="__main__":
    main()
