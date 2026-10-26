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

def assert_expected_receipt(current:dict,expected:dict)->None:
    """Reject changed specimen-level MYB12 and 16 transcript-source IDs.

    The durable frozen summary is independently assembled from the original
    three public source files. Checking an expected file does not substitute
    for downloading and hashing original OSF bytes on each live run.
    """
    obs=expected["our_descriptive_audit"]
    expected_ids=[r["id"] for r in obs["strict_perfect_source_entries"]]
    checks=[
        current["n_samples_blue"]==expected["original_n_blue"],
        current["n_samples_red"]==expected["original_n_red"],
        current["total_filtered_transcripts"]==expected["original_n_transcript_entries"],
        current["perfect_blue_over_red"]==obs["n_perfect_BLUE_GT_RED_transcript_entries"],
        current["perfect_red_over_blue"]==obs["n_perfect_RED_GT_BLUE_transcript_entries"],
        current["perfect_blue_transcript_IDs"]==expected_ids,
        current["perfect_red_transcript_IDs"]==[],
        current["selected_gene"]["pairwise_blue_over_red_AUC"]==obs["blue_over_red_pairwise_AUC"],
    ]
    for current_key,frozen_key in (
        ("blue_mean_tpm","mean_blue_TPM"),
        ("red_mean_tpm","mean_red_TPM"),
        ("blue_min_tpm","min_blue_TPM"),
        ("red_max_tpm","max_red_TPM"),
    ):
        checks.append(math.isclose(
            current["selected_gene"][current_key],obs[frozen_key],
            rel_tol=0,abs_tol=1e-12))
    expected_samples=[
        {"source_original_ID":s["sample_id"],
         "paper_phenotype":"blue" if s["source_flower_colour"]=="B" else "red",
         "TPM":s["TPM"]}
        for s in obs["original_Myb12_TPM_by_sample"]
    ]
    checks.append(current["selected_gene"]["individuals"]==expected_samples)
    if len(expected_ids)!=16 or not all(checks):
        raise ValueError("Frozen original 24-sample and 16-transcript receipt mismatch")



def extract_strict_perfect_tpm_vectors(text:str, assignments:dict[str,str], transcript_ids:list[str]):
    """Second pass limited to source-verified, *already* identified IDs.

    The caller must have audited ALL 19,184 original rows beforehand.
    No test statistic chooses the 16 candidates in this function.
    """
    requested=set(transcript_ids)
    if not requested or len(requested)!=len(transcript_ids):
        raise ValueError("Expected unique nonempty pre-identified transcript IDs")
    reader=csv.reader(io.StringIO(text))
    cols=next(reader,None)
    if not cols or cols[0]!="gene_id" or len(cols)!=len(assignments)+1:
        raise ValueError("Original TPM header or specimen count mismatch")
    ids=cols[1:]
    if len(set(ids))!=len(ids) or set(ids)!=set(assignments):
        raise ValueError("TPM specimen ID mismatch for coexpression")
    vectors={}
    for row in reader:
        if not row or row[0] not in requested:
            continue
        if len(row)!=len(cols) or row[0] in vectors:
            raise ValueError("Original transcript row shape or ID duplicate")
        try:
            values=[float(v) for v in row[1:]]
        except ValueError as exc:
            raise ValueError("Non-numeric selected TPM") from exc
        if any(not math.isfinite(v) or v<0 for v in values):
            raise ValueError("Selected transcript has negative/invalid TPM")
        vectors[row[0]]=values
    if set(vectors)!=requested:
        raise ValueError("Incomplete original strict-perfect transcript inventory")
    return ids,vectors


def _pearson(x:list[float],y:list[float]):
    """Pearson r with explicitly undefined zero residual variance.

    No log-transform, biological replication claim or pseudocount.
    Floating point ties to machine precision are not claimed exact;
    a constant profile has a zero variance and returns None.
    """
    if len(x)!=len(y) or len(x)<2:
        raise ValueError("Correlation requires matched vectors with n>=2")
    xm=math.fsum(x)/len(x)
    ym=math.fsum(y)/len(y)
    dx=[a-xm for a in x]
    dy=[b-ym for b in y]
    xx=math.fsum(v*v for v in dx)
    yy=math.fsum(v*v for v in dy)
    if xx==0 or yy==0:
        return None
    xy=math.fsum(a*b for a,b in zip(dx,dy))
    return max(-1.,min(1.,xy/math.sqrt(xx*yy)))


def coexpression_of_selected_perfect(vectors:dict[str,list[float]],
                                     selected_gene:str,
                                     sample_ids:list[str],
                                     assignments:dict[str,str])->dict:
    """Source-internal diagnostic of pooled versus within-colour correlation.

    Each transcript is independently centered within the original 12 B and
    12 R flower-colour strata before Pearson r is computed over residuals.
    This is *not* a causal coregulation assay, pathway-function measurement,
    ancestor-state inference or held-out prediction experiment.
    """
    if (not sample_ids or len(sample_ids)!=len(set(sample_ids)) or
        set(sample_ids)!=set(assignments) or
        selected_gene not in vectors):
        raise ValueError("Coexpression sample IDs or MYB12-like candidate mismatch")
    states=[assignments[k] for k in sample_ids]
    if set(states)!={"B","R"} or min(states.count("B"),states.count("R"))<2:
        raise ValueError("Require >=2 measured individuals per colour stratum")
    if len(vectors)<2:
        raise ValueError("Require at least one separately identified transcript")
    for key,vals in vectors.items():
        if (len(vals)!=len(sample_ids) or
            any(not isinstance(v,(float,int)) or not math.isfinite(v) or v<0 for v in vals)):
            raise ValueError("Missing, invalid or negative original transcript TPM")

    blue_idx=[i for i,s in enumerate(states) if s=="B"]
    red_idx=[i for i,s in enumerate(states) if s=="R"]

    def center(x):
        centered=list(map(float,x))
        for ix in (blue_idx,red_idx):
            mean=math.fsum(centered[i] for i in ix)/len(ix)
            for i in ix:
                centered[i]-=mean
        return centered

    base=list(map(float,vectors[selected_gene]))
    base_centered=center(base)
    out=[]
    for transcript in sorted(set(vectors)-{selected_gene}):
        vals=list(map(float,vectors[transcript]))
        if min(vals[i] for i in blue_idx)<=max(vals[i] for i in red_idx):
            raise ValueError("Selected comparator is not a strict blue-over-red separator")
        pooled=_pearson(base,vals)
        adjusted=_pearson(base_centered,center(vals))
        out.append({
            "transcript_id":transcript,
            "pooled_r":pooled,
            "within_color_centered_r":adjusted,
            "correlation_unit":"original backcross individual",
        })
    return {
        "reference_transcript":selected_gene,
        "original_individuals":len(sample_ids),
        "blue_individuals":len(blue_idx),
        "red_individuals":len(red_idx),
        "other_blue_perfect_transcripts":len(out),
        "pooled_positive":sum(v["pooled_r"] is not None and v["pooled_r"]>0 for v in out),
        "pooled_negative":sum(v["pooled_r"] is not None and v["pooled_r"]<0 for v in out),
        "pooled_undefined":sum(v["pooled_r"] is None for v in out),
        "color_centered_positive":sum(v["within_color_centered_r"] is not None and v["within_color_centered_r"]>0 for v in out),
        "color_centered_negative":sum(v["within_color_centered_r"] is not None and v["within_color_centered_r"]<0 for v in out),
        "color_centered_undefined":sum(v["within_color_centered_r"] is None for v in out),
        "per_transcript":out,
        "descriptive_in_sample_only":True,
        "conditioning_on_flower_color_not_molecular_intervention":True,
        "not_independent_causal_coregulation":True,
        "no_natural_regain_or_fitness_assay":True,
    }




def colour_conditioned_sensitivity(vectors:dict[str,list[float]],
                                    selected_gene:str,
                                    sample_ids:list[str],
                                    assignments:dict[str,str])->dict:
    """Leave-one-original-individual-out and log1p audit of stratified r.

    Sensitivity analysis only: the 24 backcross individuals and original
    group-selected transcripts are not independent genetic causes or
    prospective ancestral regains. Each omission re-estimates the blue/red
    means. Zero-variance group-specific correlations remain undefined.
    """
    primary=coexpression_of_selected_perfect(
        vectors,selected_gene,sample_ids,assignments)
    n=len(sample_ids)
    counts={c:list(assignments[k] for k in sample_ids).count(c) for c in ("B","R")}
    if min(counts.values())<3:
        raise ValueError("Jackknife requires at least three original individuals per colour group")
    transformations={k:[math.log1p(v) for v in value] for k,value in vectors.items()}
    logged=coexpression_of_selected_perfect(
        transformations,selected_gene,sample_ids,assignments)
    logged_lookup={z["transcript_id"]:z for z in logged["per_transcript"]}
    delete_lookups=[]
    for dropped in range(n):
        reduced_ids=[key for i,key in enumerate(sample_ids) if i!=dropped]
        reduced_vectors={k:[val for i,val in enumerate(values) if i!=dropped]
                         for k,values in vectors.items()}
        diagnostic=coexpression_of_selected_perfect(
            reduced_vectors,selected_gene,reduced_ids,
            {k:assignments[k] for k in reduced_ids})
        delete_lookups.append({z["transcript_id"]:z for z in diagnostic["per_transcript"]})

    b_indices=[i for i,key in enumerate(sample_ids) if assignments[key]=="B"]
    r_indices=[i for i,key in enumerate(sample_ids) if assignments[key]=="R"]
    reference=vectors[selected_gene]
    rows=[]
    for entry in primary["per_transcript"]:
        gene=entry["transcript_id"]
        initial=entry["within_color_centered_r"]
        transformed=logged_lookup[gene]["within_color_centered_r"]
        all_loo=[d[gene]["within_color_centered_r"] for d in delete_lookups]
        same=sum(v is not None and initial is not None and
                 (v>0)==(initial>0) and v!=0 and initial!=0
                 for v in all_loo)
        opposite=sum(v is not None and initial is not None and
                     (v>0)!=(initial>0) and v!=0 and initial!=0
                     for v in all_loo)
        undefined=sum(v is None or v==0 for v in all_loo)
        blue_r=_pearson([reference[i] for i in b_indices],
                        [vectors[gene][i] for i in b_indices])
        red_r=_pearson([reference[i] for i in r_indices],
                       [vectors[gene][i] for i in r_indices])
        rows.append({
            "transcript_id":gene,
            "original_centered_r":initial,
            "loocv_n":n,
            "omissions_same_nonzero_sign":same,
            "omissions_opposite_nonzero_sign":opposite,
            "omissions_undefined_or_zero":undefined,
            "loocv_min_r":min((v for v in all_loo if v is not None),default=None),
            "loocv_max_r":max((v for v in all_loo if v is not None),default=None),
            "log1p_centered_r":transformed,
            "log1p_same_nonzero_sign":(
                None if initial is None or transformed is None or
                        initial==0 or transformed==0
                else (initial>0)==(transformed>0)),
            "blue_only_r":blue_r,
            "red_only_r":red_r,
        })
        if same+opposite+undefined!=n:
            raise ArithmeticError("Leave-one-out sign partition lost specimens")
    neg=[row for row in rows if row["original_centered_r"] is not None
         and row["original_centered_r"]<0]
    return {
        "original_individuals":n,
        "original_blue":counts["B"],
        "original_red":counts["R"],
        "other_preselected_transcript_entries":len(rows),
        "raw_color_centered_negative":len(neg),
        "negative_with_all_omissions_same_sign":sum(
            x["omissions_same_nonzero_sign"]==n for x in neg),
        "negative_with_at_least_one_omission_sign_reversal":sum(
            x["omissions_opposite_nonzero_sign"]>0 for x in neg),
        "negative_preserved_under_log1p":sum(
            x["log1p_same_nonzero_sign"] is True for x in neg),
        "original_groups_recentered_after_each_omission":True,
        "not_independent_holdout_or_predictive_cross_validation":True,
        "not_a_formal_uncertainty_interval_or_causal_regulatory_assay":True,
        "one_original_backcross_no_ancestral_loss_regain_or_fitness":True,
        "per_transcript":rows,
    }


def assert_expected_coexpression(current:dict,expected_receipt:dict)->None:
    """Fail closed on drift from the independently frozen source-derived result.

    The frozen values were produced by a previous MD5-gated OSF live replay,
    and are descriptive within a single segregating original backcross.
    """
    if (expected_receipt.get("source_original_tpm_md5")!=SOURCE["tpm"]["md5"] or
        expected_receipt.get("source_original_manifest_md5")!=SOURCE["sample"]["md5"] or
        expected_receipt.get("source_original_gene_crosswalk_md5")!=SOURCE["gene_crosswalk"]["md5"]):
        raise ValueError("Frozen original source digest mismatch")
    expected=expected_receipt.get("source_verified_descriptive_result",{})
    fields=("reference_transcript","original_individuals","blue_individuals",
            "red_individuals","other_blue_perfect_transcripts",
            "pooled_positive","pooled_negative","pooled_undefined",
            "color_centered_positive","color_centered_negative",
            "color_centered_undefined","descriptive_in_sample_only",
            "conditioning_on_flower_color_not_molecular_intervention",
            "not_independent_causal_coregulation",
            "no_natural_regain_or_fitness_assay")
    if any(current.get(key)!=expected.get(key) for key in fields):
        raise ValueError("Frozen phenotype-conditioned coexpression counts mismatch")
    observed=current.get("per_transcript",[])
    frozen=expected.get("per_transcript",[])
    if len(observed)!=15 or len(frozen)!=15:
        raise ValueError("Frozen coexpression must contain 15 original comparators")
    for a,b in zip(observed,frozen):
        if (a.get("transcript_id")!=b.get("transcript_id") or
            a.get("correlation_unit")!=b.get("correlation_unit")):
            raise ValueError("Frozen transcript identity or measurement unit mismatch")
        for name in ("pooled_r","within_color_centered_r"):
            x,y=a.get(name),b.get(name)
            if x is None or y is None:
                if x is not None or y is not None:
                    raise ValueError("Frozen undefined correlation mismatch")
            elif not (math.isfinite(x) and math.isfinite(y) and
                      math.isclose(x,y,abs_tol=1e-12,rel_tol=0)):
                raise ValueError("Frozen per-transcript coexpression drift")


def report():
    source_bytes={name:fetch_bytes(meta) for name,meta in SOURCE.items()}
    originals={name:val.decode("utf-8-sig") for name,val in source_bytes.items()}
    manifest=parse_manifest(originals["sample"])
    x=audit_tpm_table(originals["tpm"],manifest)
    x["author_DE_source"]=crosscheck_author_DE(originals["gene_crosswalk"])
    selected_ids=x["perfect_blue_transcript_IDs"]
    ids,original_perfect=extract_strict_perfect_tpm_vectors(
        originals["tpm"],manifest,selected_ids)
    x["color_conditioned_coexpression"]=coexpression_of_selected_perfect(
        original_perfect,MYB_ID,ids,manifest)
    x["color_conditioned_sensitivity"]=colour_conditioned_sensitivity(
        original_perfect,MYB_ID,ids,manifest)
    if x["color_conditioned_coexpression"]["other_blue_perfect_transcripts"]!=15:
        raise ValueError("Coexpression source-selected comparator count is not 15")
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
            "F3primeH_Myb12_0p91_original_correlation_not_independently_replicated_from_filtered_tpm_matrix":True,
            "Cirsium_young_leaf_source_is_not_a_flower_pigment_comparison":True,
        }
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--expected",type=Path)
    parser.add_argument("--conditioned-expected",type=Path)
    args=parser.parse_args()
    result=report()
    if args.expected:
        frozen=json.loads(args.expected.read_text(encoding="utf-8"))
        assert_expected_receipt(result,frozen)
    if args.conditioned_expected:
        audit=json.loads(args.conditioned_expected.read_text(encoding="utf-8"))
        assert_expected_coexpression(result["color_conditioned_coexpression"],audit)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("IOCHROMA_COLOR_SENSITIVITY_DETAIL",json.dumps(
        result["color_conditioned_sensitivity"],sort_keys=True))
    print("IOCHROMA_COLOR_CONDITIONED_DETAIL",json.dumps(
        result["color_conditioned_coexpression"],sort_keys=True))
    print("IOCHROMA_REAL_SOURCE_ADMISSION_PASS",json.dumps({
        "n_samples":result["n_samples_blue"]+result["n_samples_red"],
        "transcripts":result["total_filtered_transcripts"],
        "perfect_blue":result["perfect_blue_over_red"],
        "perfect_red":result["perfect_red_over_blue"],
        "MYB12_blue_mean":result["selected_gene"]["blue_mean_tpm"],
        "MYB12_red_mean":result["selected_gene"]["red_mean_tpm"],
        "coexpression_pooled_positive":result["color_conditioned_coexpression"]["pooled_positive"],
        "coexpression_within_color_negative":result["color_conditioned_coexpression"]["color_centered_negative"],
        "coexpression_within_color_undefined":result["color_conditioned_coexpression"]["color_centered_undefined"],
        "claim_limit":"Original backcross expression reanalysis, not historical evolution",
    },sort_keys=True))

if __name__=="__main__":main()
