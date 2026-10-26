"""Reanalyse REAL 2025 Iochroma red/blue floral RNA expression at source grain.

Original authors: Wheeler, Larter & Smith (2025), G3, DOI
10.1093/g3journal/jkaf230; OSF https://osf.io/j5m8f/ .
This workflow only recovers and re-audits the authors' 24 segregating
backcross individuals (12 red, 12 blue); it does not identify genetic
causality, ancestral trait loss, evolutionary re-expression or fixation.

Source files are re-downloaded from official OSF IDs and must match the
original provider's exact source MD5 and byte lengths. No raw 4.8-MB
original transcriptome is rehosted in this repository.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import sys
import time
import urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OSF_ROOT="https://osf.io/j5m8f/"
PAPER_DOI="10.1093/g3journal/jkaf230"
SOURCE={
    "tpm":{"id":"6622e385ea2b630b23532fe5","bytes":4792342,
           "md5":"dcf46c27169d1b84612a2fd4c43d4642",
           "path":"WGCNA/transcript-by-species-filtered-TPM-matrix-transpose.csv"},
    "sample":{"id":"6622e067387f9d08bf7ac834","bytes":955,
              "md5":"a9264d3a09fca80efe62ce3eaee32046",
              "path":"differental_expression_analysis/samples.csv"},
    "gene_crosswalk":{"id":"6622e4e220a88e09c136efe0","bytes":6780,
              "md5":"d6584d58743b0b5f901ae0fe030488d4",
              "path":"differental_expression_analysis/de-transcripts-coordinates-annotations.csv"}
}
MYB_ID="TRINITY_DN3000_c0_g1_i10.p1"
EXPECTED_LOSS_STUDY="10.1002/dvdy.82"

def validate_payload(raw:bytes,metadata:dict)->None:
    if len(raw)!=metadata["bytes"] or hashlib.md5(raw).hexdigest()!=metadata["md5"]:
        raise ValueError("OSF original source bytes or MD5 changed: refuse without re-audit")

def fetch_bytes(source:dict,retries=2)->bytes:
    url="https://osf.io/download/"+source["id"]+"/"
    for attempt in range(retries):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"ANTO-original-Iochroma-floral-source/1.0",
                "Accept":"text/csv,text/plain;q=0.9,*/*;q=0.8"})
            with urllib.request.urlopen(req,timeout=55) as response:
                buf=response.read(source["bytes"]+1)
            validate_payload(buf,source)
            return buf
        except Exception:
            if attempt+1>=retries:
                raise
            time.sleep(1)
    raise RuntimeError("OSF acquisition failed")

def parse_manifest(text:str)->dict[str,str]:
    reader=csv.DictReader(io.StringIO(text))
    if not reader.fieldnames or not {"ID","condition","path","location"}.issubset(reader.fieldnames):
        raise ValueError("Original OSF sample manifest fields missing")
    assignments={}
    for row in reader:
        id_=row["ID"]
        if not id_ or id_ in assignments:
            raise ValueError("Original specimen ID missing or duplicate")
        state=row["condition"]
        if state not in ("B","R"):
            raise ValueError("Unexpected flower color source state")
        assignments[id_]=state
    if not assignments or not set(assignments.values())=={"B","R"}:
        raise ValueError("Original specimen color contrast incomplete")
    return assignments

def audit_tpm_table(text:str,assignments:dict[str,str],selected_gene=MYB_ID)->dict:
    reader=csv.reader(io.StringIO(text))
    try:
        columns=next(reader)
    except StopIteration:
        raise ValueError("Original source TPM matrix empty")
    if not columns or columns[0]!="gene_id":
        raise ValueError("Expected original gene_id first field")
    sample_ids=columns[1:]
    if (len(sample_ids)!=len(set(sample_ids)) or
        set(sample_ids)!=set(assignments) or not sample_ids):
        raise ValueError("TPM sample IDs do not match original 24 specimen sample labels")
    group_B=[i for i,k in enumerate(sample_ids) if assignments[k]=="B"]
    group_R=[i for i,k in enumerate(sample_ids) if assignments[k]=="R"]
    if not group_B or not group_R:
        raise ValueError("Missing B/R sample group")
    ids=set()
    n=perfect_blue=perfect_red=0
    perfect_ids_blue=[]
    perfect_ids_red=[]
    selected=None
    for row in reader:
        if not row:
            continue
        id_=row[0]
        if not id_ or id_ in ids:
            raise ValueError("gene_id duplicate or empty in original transcript matrix")
        ids.add(id_)
        if len(row)!=len(columns):
            raise ValueError("Source TPM row has wrong sample count")
        try:
            vals=list(map(float,row[1:]))
        except ValueError as e:
            raise ValueError("Invalid non-numeric TPM") from e
        if any(not math.isfinite(v) or v<0 for v in vals):
            raise ValueError("Negative, missing, NaN, or infinite TPM")
        blue=[vals[i] for i in group_B]
        red=[vals[i] for i in group_R]
        sepB=min(blue)>max(red)
        sepR=min(red)>max(blue)
        n+=1
        perfect_blue+=int(sepB)
        perfect_red+=int(sepR)
        if sepB: perfect_ids_blue.append(id_)
        if sepR: perfect_ids_red.append(id_)
        if id_==selected_gene:
            wins=sum(1. if b>r else .5 if b==r else 0.
                     for b in blue for r in red)
            selected={
                "original_transcript_id":id_,
                "blue_mean_tpm":sum(blue)/len(blue),
                "red_mean_tpm":sum(red)/len(red),
                "blue_min_tpm":min(blue),
                "red_max_tpm":max(red),
                "blue_max_tpm":max(blue),
                "red_min_tpm":min(red),
                "strict_rank_discrimination":sepB,
                "pairwise_blue_over_red_AUC":wins/(len(blue)*len(red)),
                "individuals":[{
                    "source_original_ID":k,
                    "paper_phenotype":"blue" if assignments[k]=="B" else "red",
                    "TPM":vals[i]
                } for i,k in enumerate(sample_ids)],
                "in_sample_and_author_prior_selected":True,
                "no_unique_genetic_causality":True
            }
    if selected is None:
        raise ValueError(f"Selected MYB12-like original transcript {selected_gene} not found")
    return {
        "total_filtered_transcripts":n,
        "n_samples_blue":len(group_B),
        "n_samples_red":len(group_R),
        "perfect_blue_over_red":perfect_blue,
        "perfect_red_over_blue":perfect_red,
        "perfect_transcript_fraction":(perfect_blue+perfect_red)/n,
        "perfect_blue_transcript_IDs":sorted(perfect_ids_blue),
        "perfect_red_transcript_IDs":sorted(perfect_ids_red),
        "selected_gene":selected,
        "trait_loss_then_historical_regain_measured":False,
        "causal_MYB_gene_demonstrated_by_reanalysis":False,
        "white_morph_vs_fixed_white_species_measured":False,
        "single_original_backcross_not_independent_macro_clades":True,
    }

def crosscheck_author_DE(text:str)->dict:
    reader=csv.DictReader(io.StringIO(text))
    required={"Transcript","Icyan_chrom","log2FoldChange","padj"}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise ValueError("Missing original gene candidate provenance fields")
    rows=[x for x in reader if x["Transcript"]==MYB_ID]
    if len(rows)!=1:
        raise ValueError("Unexpected source candidate transcript ambiguity")
    row=rows[0]
    lfc=float(row["log2FoldChange"])
    if not math.isclose(lfc,-6.35096062795051,rel_tol=0,abs_tol=1e-10):
        raise ValueError("Original author MYB12-like DE source annotation drift")
    if row["Icyan_chrom"]!="IC1.0ch05":
        raise ValueError("Original MYB12 chromosome annotation drift")
    return {
        "transcript":MYB_ID,"chromosome":row["Icyan_chrom"],
        "source_log2FoldChange_red_relative_blue":lfc,
        "original_source_padj":float(row["padj"]),
        "source_note":"Authors' pre-existing DESeq2 effect; cannot treat TPM mean ratio as an independently fitted DESeq2 coefficient",
    }

def report():
    source_bytes={name:fetch_bytes(meta) for name,meta in SOURCE.items()}
    originals={name:val.decode("utf-8-sig") for name,val in source_bytes.items()}
    manifest=parse_manifest(originals["sample"])
    x=audit_tpm_table(originals["tpm"],manifest)
    x["author_DE_source"]=crosscheck_author_DE(originals["gene_crosswalk"])
    if (x["total_filtered_transcripts"]!=19184 or x["n_samples_red"]!=12 or
        x["n_samples_blue"]!=12 or x["perfect_blue_over_red"]!=16 or
        x["perfect_red_over_blue"]!=0):
        raise ValueError("Original author source population or audit denominator drift")
    chosen=x["selected_gene"]
    if (not math.isclose(chosen["blue_mean_tpm"],11.824398166666667,abs_tol=1e-9) or
        not math.isclose(chosen["red_mean_tpm"],.16201,abs_tol=1e-9) or
        chosen["pairwise_blue_over_red_AUC"]!=1.0):
        raise ValueError("Source original MYB12 line and phenotype contrast changed")
    return {
        "schema":"anto_iochroma_2025_original_backcross_myb12_flower_TPM_audit_v1",
        "audit_date":"2026-10-26",
        "source":"Wheeler et al 2025 original OSF j5m8f",
        "original_publication_doi":PAPER_DOI,
        "dataset_url":OSF_ROOT,
        "source_files":{
            name:{"osf_file_id":meta["id"],"registered_path":meta["path"],
                  "expected_md5":meta["md5"],
                  "expected_bytes":meta["bytes"],
                  "original_source_hash_verified":True}
            for name,meta in SOURCE.items()
        },
        "related_fixed_vs_polymorphic_white_flower_study":{
            "doi":EXPECTED_LOSS_STUDY,
            "actual_2019_raw_individual_data_reanalyzed":False,
            "dataset_DOI":"10.5061/dryad.p5dq84v",
            "current_raw_archive_access":"HTTP_401_FOR_2019_DRYAD_API",
            "original_paper_claim":"white morph and fixed white species have distinct gene regulation despite similarly absent flower anthocyanins",
        },
        **x,
        "claim_boundaries":{
            "24_original_specimens_are_one_backcross_not_24_independent_evolutionary_origins":True,
            "19184_transcript_rows_are_not_19184_independent_genes_or_tests":True,
            "16_in_sample_perfect_separators_are_not_16_causal_loci":True,
            "MYB12_like_candidates_were_selected_by_original_authors_before_this_audit":True,
            "AUC_1_in_backcross_is_not_prospective_clade_holdout_accuracy":True,
            "no_loss_ancestry_phenotype_regain_or_seed_fitness_measured":True,
            "F3H_Myb12_0p91_original_correlation_not_independently_replicated_from_filtered_tpm_matrix":True,
            "Cirsium_young_leaf_source_is_not_a_flower_pigment_comparison":True,
        }
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result=report()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("IOCHROMA_REAL_SOURCE_ADMISSION_PASS",json.dumps({
        "n_samples":result["n_samples_blue"]+result["n_samples_red"],
        "transcripts":result["total_filtered_transcripts"],
        "perfect_blue":result["perfect_blue_over_red"],
        "perfect_red":result["perfect_red_over_blue"],
        "MYB12_blue_mean":result["selected_gene"]["blue_mean_tpm"],
        "MYB12_red_mean":result["selected_gene"]["red_mean_tpm"],
        "claim_limit":"Original backcross expression reanalysis, not historical evolution",
    },sort_keys=True))

if __name__=="__main__":main()
