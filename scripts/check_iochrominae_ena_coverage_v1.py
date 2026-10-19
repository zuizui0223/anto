"""Reproducible post-retrieval coverage audit of 58 original Iochrominae records.

Reads the actual source-only ENA receipt, derives per-species deposited
accession coverage by *registered accession block*, validates independently
retained publication-scope ledger, and stops if original provenance drifts.

None of the accession blocks is a verified functional gene deletion,
evolutionary reacquisition, enzyme activity or trait-fitness observation.
"""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

SOURCE_RECEIPT_SHA256 = "f528ca11e0ded70049d3def146b03b86058c4916a40fb35be59089052590e602"

def classify_number(acc):
    try:
        n=int(acc.removeprefix("KT"))
    except (AttributeError,ValueError) as exc:
        raise ValueError("Unexpected original accession") from exc
    if not acc.startswith("KT") or n<898394 or n>898451:
        raise ValueError("Accession outside author-registered KT898394-KT898451 span")
    return ("F3h" if n<=898415 else "Chi" if n<=898437 else "Dfr accession block")

def original_source_coverage(report):
    rows=report.get("genbank_record_provenance",[])
    if report.get("source_paper_DOI")!="10.1186/s12862-016-0675-3":
        raise ValueError("Mismatched original scientific source DOI")
    if report.get("expected_original_accession_count")!=58 or len(rows)!=58:
        raise ValueError("Original 58-record source incomplete")
    ids=[r["accession"] for r in rows]
    if len(set(ids))!=58 or set(ids)!={f"KT{i}" for i in range(898394,898452)}:
        raise ValueError("Missing/repeated accession")
    groups=defaultdict(lambda:{"F3h":0,"Chi":0,"Dfr accession block":0,
                               "Dfr accession block without CDS":False})
    length_by_group=defaultdict(list)
    products=Counter()
    missing_CDS=[]
    for row in rows:
        species=row["organism_as_submitted"]
        if not isinstance(species,str) or not species.strip():
            raise ValueError("Unadmitted biological taxon ID")
        group=classify_number(row["accession"])
        groups[species][group]+=1
        n=row["sequence_length_nt"]
        if type(n) is not int or n<500 or n>2500:
            raise ValueError("Unexpected accession sequence length")
        length_by_group[group].append(n)
        cds=row["cds_annotations"]
        if len(cds)!=row["cds_feature_count"] or len(cds)>1:
            raise ValueError("CDS annotations inconsistent")
        if not cds:
            missing_CDS.append({"accession":row["accession"],"species":species})
            if group!="Dfr accession block":
                raise ValueError("Previously annotated F3h/Chi unexpectedly lost CDS")
            groups[species]["Dfr accession block without CDS"]=True
        else:
            if not cds[0]["translation_annotated"]:
                raise ValueError("Original CDS translation annotation absent")
            product=cds[0]["product_annotations"]
            correct={
                "F3h":"flavanone 3-hydroxylase",
                "Chi":"chalcone flavanone isomerase",
                "Dfr accession block":"dihydroflavonol reductase"
            }[group]
            if product != [correct]:
                raise ValueError(f"Gene product inconsistency in {row['accession']}")
            products[correct]+=1
        if row.get("functional_enzyme_activity_measured") is not False or row.get("real_regained_flower_colour_measured") is not False:
            raise ValueError("Unexpected promotion to biological function/reversion")
        if len(row.get("record_file_sha256",""))!=64 or len(row.get("sequence_sha256",""))!=64:
            raise ValueError("Source hashes missing")

    if len(groups)!=22 or sum(p["F3h"]==1 and p["Chi"]==1 for p in groups.values())!=22:
        raise ValueError("Original 22 submitted species source missing")
    if products!={"flavanone 3-hydroxylase":22,"chalcone flavanone isomerase":22,"dihydroflavonol reductase":12}:
        raise ValueError("Product annotation coverage changed")
    if len(missing_CDS)!=2:
        raise ValueError("Missing CDS features changed")
    if sum(v["Dfr accession block"]==1 for v in groups.values())!=14:
        raise ValueError("Newly deposited Dfr accession coverage changed")
    compact_species=[
        {"submitted_species_name":sp,
         "has_new_F3h_record":True,
         "has_new_Chi_record":True,
         "has_new_Dfr_accession_in_this_span":v["Dfr accession block"]==1,
         "Dfr_accession_missing_cds_feature":v["Dfr accession block without CDS"],
         "no_Dfr_in_this_new_accession_span_is_gene_loss":False}
        for sp,v in sorted(groups.items())
    ]
    counts=Counter(r["sequence_length_nt"] for r in rows)
    if counts != Counter({1089:22,705:22,1146:11,1145:1,1143:1,1170:1}):
        raise ValueError("Original nucleotide lengths differ from recorded source")
    return {"per_taxon":compact_species,
            "summary":{
                "source_accession_count":58,
                "source_taxon_names":22,
                "F3h_accessions":len(length_by_group["F3h"]),
                "Chi_accessions":len(length_by_group["Chi"]),
                "Dfr_new_accessions":len(length_by_group["Dfr accession block"]),
                "Dfr_product_annotated":products["dihydroflavonol reductase"],
                "Dfr_new_acc_without_CDS":len(missing_CDS),
                "source_records_with_ambiguous_base":sum(bool(r["ambiguous_base_count"]) for r in rows),
                "all_ambiguous_base_positions":sum(r["ambiguous_base_count"] for r in rows),
                "full_three_locus_new_accession_species":14,
                "Dfr_new_accession_missing_species":8,
                "cds_record_count":sum(r["cds_feature_count"] for r in rows),
                "source_NCBI_record_sequence_length_min":min(counts),
                "source_NCBI_record_sequence_length_max":max(counts),
                "does_not_identify_absent_Dfr_genes_or_gene_function":True,
            },
            "no_CDS_records":sorted(missing_CDS,key=lambda x:x["accession"])}

def enforce_archive(receipt, expected):
    out=original_source_coverage(receipt)
    if out["per_taxon"] != expected["per_taxon"]:
        raise ValueError("Committed exact species record coverage differs")
    e=expected["observed"]
    if (out["summary"]["source_accession_count"]!=e["successful_accessions"]
        or out["summary"]["source_taxon_names"]!=e["submitted_taxon_names"]
        or out["summary"]["Dfr_new_accessions"]!=e["source_accession_family_groups"][2]["n"]
        or out["summary"]["full_three_locus_new_accession_species"]!=e["taxa_with_all_three_accession_groups"]
        or out["summary"]["cds_record_count"]!=e["annotated_CDS_features"]
        or out["summary"]["source_records_with_ambiguous_base"]!=e["records_with_ambiguous_nucleotides"]
        or out["summary"]["all_ambiguous_base_positions"]!=e["ambiguous_base_positions"]
        or out["no_CDS_records"]!=sorted(e["source_no_CDS_feature_accessions"],key=lambda x:x["accession"])):
        raise ValueError("Committed original-public-accession aggregate source drift")
    if expected["independent_real_reexpression_rate_estimate"] is not False:
        raise ValueError("False biological revival claim")
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--expected",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    raw=a.source.read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=SOURCE_RECEIPT_SHA256:
        raise ValueError(f"Original 58-record source receipt hash drift: {digest}")
    report=json.loads(raw)
    expected=json.loads(a.expected.read_text(encoding="utf-8"))
    out=enforce_archive(report,expected)
    out["source_receipt_sha256"]=digest
    out["source_description"]="58 original INSDC accessions and deposited CDS annotations ONLY"
    out["biological_return_functionality_inferred"]=False
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("IOCHROMINAE_58_SOURCE_COVERAGE_VERIFIED",json.dumps(out["summary"],sort_keys=True))

if __name__=="__main__":main()
