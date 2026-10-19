"""Real public GenBank/ENA accession audit for Ho & Smith (2016) coding data.

DOI 10.1186/s12862-016-0675-3 lists KT898394-KT898451.
Retrieve EMBL records via the ENA Browser API, not the Dryad 401 download.
The original authors' intact-gene and relaxed constraint conclusions are
not independently re-proven by the mere presence of sequence records.
Never infer an organism's visible phenotype, ancestral regain, floral protein
activity, matched sampled individual, dN/dS or ecological fitness from this.
Raw records are NOT rehosted; record-level SHA and annotations only.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from urllib.request import Request,urlopen

from Bio import SeqIO

ACCESSIONS=tuple(f"KT{i}" for i in range(898394,898452))
SOURCE_DOI="10.1186/s12862-016-0675-3"
BASE="https://www.ebi.ac.uk/ena/browser/api/embl/"
MAX_RECORD_BYTES=120_000


def get_record(accession,timeout=30):
    if accession not in ACCESSIONS:
        raise ValueError("Unregistered accession, stop before source expansion")
    url=BASE+accession
    req=Request(url,headers={"User-Agent":"ANTO-source-qualification/1.0",
                             "Accept":"text/plain"})
    with urlopen(req,timeout=timeout) as res:
        b=res.read(MAX_RECORD_BYTES+1)
    if len(b)>MAX_RECORD_BYTES or not b:
        raise ValueError("Missing or oversized raw GenBank/ENA record")
    return b


def parse_record(expected_accession,record_bytes):
    raw=record_bytes.decode("utf-8-sig")
    records=list(SeqIO.parse(io.StringIO(raw),"embl"))
    if len(records)!=1:
        raise ValueError("An EMBL API record must contain exactly one accession")
    record=records[0]
    current=record.id.split(".")[0]
    if current!=expected_accession or not re.fullmatch(r"KT\d{6}",current):
        raise ValueError("Record accession mismatches preregistered accession")
    dna=str(record.seq).upper()
    if not dna or set(dna)-set("ACGTNRYKMSWBDHV"):
        raise ValueError("Sequence contains invalid bases or is missing")
    cds=[]
    for feature in record.features:
        if feature.type=="CDS":
            cds.append({
                "gene_annotations":feature.qualifiers.get("gene",[]),
                "product_annotations":feature.qualifiers.get("product",[]),
                "protein_accessions":feature.qualifiers.get("protein_id",[]),
                "translation_annotated":bool(feature.qualifiers.get("translation")),
                "partial_location_tag":str(feature.location).find("<")>=0 or str(feature.location).find(">")>=0,
            })
    return {
        "accession":expected_accession,
        "record_id":record.id,
        "record_version":record.id[len(current):].lstrip("."),
        "organism_as_submitted":record.annotations.get("organism"),
        "organism_recorded_not_phenotype_verified":True,
        "sequence_length_nt":len(dna),
        "ambiguous_base_count":sum(c not in "ACGT" for c in dna),
        "sequence_sha256":hashlib.sha256(dna.encode("ascii")).hexdigest(),
        "record_file_sha256":hashlib.sha256(record_bytes).hexdigest(),
        "cds_feature_count":len(cds),
        "cds_annotations":cds,
        "functional_enzyme_activity_measured":False,
        "real_regained_flower_colour_measured":False,
    }


def retrieve_all():
    results=[]
    errors=[]
    def process(acc):
        return parse_record(acc,get_record(acc))
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(process,acc):acc for acc in ACCESSIONS}
        for future in as_completed(futures):
            acc=futures[future]
            try:
                results.append(future.result())
            except Exception as exc:
                errors.append({"accession":acc,
                               "source_failure_type":type(exc).__name__,
                               "source_failure_message":str(exc)[:180]})
    results.sort(key=lambda a:a["accession"])
    errors.sort(key=lambda a:a["accession"])
    return results,errors


def summarize(results,errors):
    if {r["accession"] for r in results}&{e["accession"] for e in errors}:
        raise ValueError("Duplicate and failed accession overlap")
    if len(results)+len(errors)!=len(ACCESSIONS):
        raise ValueError("Missing accession from admission denominator")
    annotated_genes={}
    for record in results:
        for cds in record["cds_annotations"]:
            for label in cds["gene_annotations"]:
                annotated_genes[label]=annotated_genes.get(label,0)+1
    return {
        "schema":"anto_iochrominae_published_genbank_accession_audit_v1",
        "source_paper_DOI":SOURCE_DOI,
        "original_accessions_span":"KT898394-KT898451",
        "source_repository":"INSDC/European Nucleotide Archive",
        "api_endpoint":"https://www.ebi.ac.uk/ena/browser/api/embl/{accession}",
        "scientific_admission":"SOURCE_RECORD_AVAILABILITY_ONLY",
        "expected_original_accession_count":58,
        "retrieved_count":len(results),
        "not_retrieved_count":len(errors),
        "all_records_retrieved":not errors,
        "genbank_record_provenance":results,
        "fail_closed_accessions":errors,
        "distinct_submitted_organism_names":len(set(x["organism_as_submitted"] for x in results if x["organism_as_submitted"])),
        "reported_CDS_gene_annotations":annotated_genes,
        "original_author_coding_integrity_result_not_independently_retested":True,
        "observed_flower_colour_state_from_2016_source_joined":False,
        "published_historical_regain_event_count_reattributed":False,
        "molecular_enzyme_activity_assayed_by_this_project":False,
        "flower_trait_regain_and_persistence_data_collected":False,
        "raw_sequences_redistributed":False,
        "no_false_independent_clade_or_transition_replication":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    results,errors=retrieve_all()
    report=summarize(results,errors)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("IOCHROMINAE_INSDC_SOURCE_RETRIEVAL",json.dumps({
        "expected":58,"retrieved":len(results),"unavailable":len(errors),
        "distinct_submitted_taxa":report["distinct_submitted_organism_names"],
        "gene_annotation_counts":report["reported_CDS_gene_annotations"],
        "not_a_trait_regain_measurement":True
    },sort_keys=True))


if __name__=="__main__":
    main()
