"""Source-linked original GenBank accession audit for Ho & Smith 2016.

Original article: Ho WW & Smith SD, BMC Evol Biol 16:98 (2016).
Article's availability section registers KT898394–KT898451 (58 accessions).
This replay uses NCBI Entrez/GenBank registered accession records, NOT
the previously restricted Dryad archive. No original sequences are rehosted.

The primary outcome is an independent machine record of which exact published
accessions are retrievable, their source taxon/gene metadata and explicit
coding annotations. This is NOT an analysis of causal flower regain, and
retained coding annotation is not a functional enzymatic validation.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
import time
import urllib.request
from pathlib import Path
from Bio import SeqIO

FIRST=898394
LAST=898451
ACCESSIONS=tuple(f"KT{i}" for i in range(FIRST,LAST+1))
SOURCE_DOI="10.1186/s12862-016-0675-3"
ARTICLE_URL="https://link.springer.com/article/"+SOURCE_DOI
API="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def accession_query(batch):
    if not batch or any(x not in ACCESSIONS for x in batch):
        raise ValueError("Only exact author-registered KT898394–KT898451 accessions")
    return (API+"?db=nuccore&rettype=gb&retmode=text&tool=ANTO-iochrominae-2016"
            +"&id="+",".join(batch))

def fetch_batch(batch, opener=None):
    url=accession_query(batch)
    req=urllib.request.Request(url,headers={"User-Agent":"ANTO-original-accessions/1.0 (publication provenance audit)"})
    if opener is None: opener=urllib.request.urlopen
    with opener(req,timeout=60) as response:
        raw=response.read(1_600_001)
    if not raw or len(raw)>1_600_000:
        raise RuntimeError("NCBI_EMPTY_OR_UNEXPECTED_BATCH_SIZE")
    try:
        records=list(SeqIO.parse(io.StringIO(raw.decode("utf-8")),"genbank"))
    except (UnicodeError,ValueError) as ex:
        raise RuntimeError("NCBI_NOT_PARSEABLE_GENBANK") from ex
    obs=[r.id.split(".")[0] for r in records]
    if len(obs)!=len(batch) or set(obs)!=set(batch) or len(obs)!=len(set(obs)):
        raise RuntimeError("NCBI_ACCESSION_RECORD_ID_MISMATCH")
    return records

def gene_label(record):
    text=" ".join([record.description,*(
        q for feat in record.features
        for field in ("gene","product","note")
        for q in feat.qualifiers.get(field,[])
    )]).lower()
    patterns={
        "CHI":r"\bchi\b|chalcone\s+isomerase",
        "F3H":r"\bf3h\b|flavanone[ -]?3[- ]?hydroxylase|flavanone[ -]?3[- ]?dioxygenase",
        "DFR":r"\bdfr\b|dihydroflavonol[ -]?(?:4[- ]?)?reductase",
    }
    names=[g for g,p in patterns.items() if re.search(p,text)]
    return names[0] if len(names)==1 else ("MULTIPLE_GENE_ALIASES" if names else "UNRESOLVED_GENE_LABEL")

def compact_record(rec):
    accession=rec.id.split(".")[0]
    if accession not in ACCESSIONS:
        raise ValueError("Unexpected accession in original author source")
    organism=rec.annotations.get("organism","")
    if not organism or organism.lower().startswith("unknown"):
        raise RuntimeError("NCBI_ORGANISM_MISSING")
    cds=[f for f in rec.features if f.type=="CDS"]
    accession_versions=[rec.id]
    # Annotation only; do NOT infer protein function from an uninterrupted ORF.
    metadata={
        "accession":accession,
        "source_record_version":rec.id,
        "organism":organism,
        "named_gene":gene_label(rec),
        "source_nucleotide_bases":len(rec.seq),
        "source_nucleotide_sequence_sha256":hashlib.sha256(str(rec.seq).upper().encode()).hexdigest(),
        "CDS_annotations":len(cds),
        "CDS_qualifier_gene":[f.qualifiers.get("gene",[]) for f in cds],
        "CDS_qualifier_product":[f.qualifiers.get("product",[]) for f in cds],
        "full_functional_anthocyanin_pathway_verified":False,
        "phenotypic_state_from_this_record":"UNKNOWN",
        "trait_reversal_event_from_this_record":"NOT_IDENTIFIABLE"
    }
    return metadata

def audit(records):
    obs=[r.id.split(".")[0] for r in records]
    if len(records)!=len(ACCESSIONS) or len(set(obs))!=len(ACCESSIONS) or set(obs)!=set(ACCESSIONS):
        raise ValueError("The original accession range is not exhaustively retrievable")
    items=sorted([compact_record(r) for r in records],key=lambda x:int(x["accession"][2:]))
    taxa=sorted(set(x["organism"] for x in items))
    genes={}
    for r in items:genes[r["named_gene"]]=genes.get(r["named_gene"],0)+1
    return {
        "schema":"anto_iochrominae_ho_smith_2016_genbank_registered_accessions_v1",
        "original_authors":"Winnie W. Ho and Stacey D. Smith",
        "original_publication_DOI":SOURCE_DOI,
        "accession_first":ACCESSIONS[0],
        "accession_last":ACCESSIONS[-1],
        "registered_accession_count":len(ACCESSIONS),
        "source_records_recovered":len(items),
        "unique_organism_names_in_deposited_records":len(taxa),
        "taxa_named_in_original_2016_paper_as_ingroup":22,
        "taxa_named_in_original_2016_paper_as_outgroups":6,
        "source_record_gene_label_counts":genes,
        "organism_names_in_retrieved_records":taxa,
        "taxon_gene_accession_inventory":{
            organism:{gene:sorted([r["accession"] for r in items
                        if r["organism"]==organism and r["named_gene"]==gene])
                     for gene in ("CHI","F3H","DFR")}
            for organism in taxa
        },
        "original_sequences_embedded_in_result":False,
        "coding_integrity_verified_by_functional_experiment":False,
        "phenotypic_color_crosswalk_independently_qualified":False,
        "phylogenetic_tree_or_dnds_originally_recomputed":False,
        "regain_or_perpetuation_probability_estimated":False,
        "not_independent_evolutionary_trait_reversal_events":True,
        "accessions":items,
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",required=True,type=Path)
    parser.add_argument("--batch-size",type=int,default=10)
    args=parser.parse_args()
    if not 1<=args.batch_size<=15:raise ValueError("Batch size 1..15")
    all_records=[]
    for i in range(0,len(ACCESSIONS),args.batch_size):
        chunk=ACCESSIONS[i:i+args.batch_size]
        all_records.extend(fetch_batch(chunk))
        if i+args.batch_size<len(ACCESSIONS):
            time.sleep(.5)
    result=audit(all_records)
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("IOCHROMINAE_ORIGINAL_KT_ACCESSION_AUDIT",json.dumps({
        "n_accessions":result["source_records_recovered"],
        "n_taxa":result["unique_organism_names_in_deposited_records"],
        "gene_labels":result["source_record_gene_label_counts"],
        "no_biological_reversion_estimated":True,
    },sort_keys=True))
    print("ORIGINAL_KT_TAXON_GENE_INVENTORY",
          json.dumps(result["taxon_gene_accession_inventory"],sort_keys=True))

if __name__=="__main__":main()
