"""Original-source biological trait/genotype crosswalk for Iochrominae 2016.

Per-species original figure 1 anthocyanin status (7 unpigmented / 15
pigmented) joined to genuine NCBI 58-accession gene *deposit* coverage,
after separately verified taxonomic synonyms.

Critically:
- GenBank CDS annotation does NOT verify enzyme function;
- No accession in KT range does NOT mean biological gene loss;
- One inferred pigment regain in I. edule is published PHYLOGENETIC
  inference, not our molecular recreation or heritable experiment;
- Two entries without a CDS feature are UNKNOWN, not loss-of-function.
"""
from pathlib import Path
import json
from collections import Counter

ROOT=Path(__file__).resolve().parents[1]
GENES=ROOT/"data/iochrominae_ho_smith_original_KT_taxon_gene_crosswalk_20261019.json"
TRAITS=ROOT/"data/iochrominae_2016_fig1_source_trait_taxon_crosswalk_v1.json"
ANNOTATIONS=ROOT/"data/iochrominae_ho_smith_2016_original_genbank_CDS_annotation_audit_20261019.json"

def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))

def calculate(gene,trait,annotation):
    matrix=gene["table"]
    species=trait["source_taxa"]
    expected=set(matrix)
    observed=[r["current_genbank_taxon"] for r in species]
    if len(species)!=22 or len(set(observed))!=22 or set(observed)!=expected:
        raise ValueError("Original 22 species or taxon concept crosswalk changed")
    status=Counter(r["floral_anthocyanin_status_from_2016_figure"] for r in species)
    if status!={"PIGMENTED":15,"UNPIGMENTED":7}:
        raise ValueError("Source published figure 1 biological status counts changed")
    missingCDS={r["accession"] for r in annotation["source_records_without_CDS_annotation"]}
    if missingCDS!={"KT898443","KT898446"}:
        raise ValueError("Source record CDS-annotation missingness has changed")
    table=[]
    for row in species:
        taxon=row["current_genbank_taxon"]
        links=matrix[taxon]
        if set(links)!={"CHI","F3H","DFR"}:
            raise ValueError("Unexpected gene key")
        for gene_name in ("CHI","F3H"):
            if len(links[gene_name])!=1:
                raise ValueError("Original CHI/F3H deposit coverage invalid")
        if len(links["DFR"])>1:
            raise ValueError("DFR multiple accession cannot be silently treated as one locus")
        dfr=links["DFR"]
        table.append({
            "current_taxon_name":taxon,
            "original_2016_figure_taxon":row["figure_2016_taxon"],
            "trait_state_original_fig1":row["floral_anthocyanin_status_from_2016_figure"],
            "taxonomic_join_type":row["nomenclature_match"],
            "CHI_accession":links["CHI"][0],
            "F3H_accession":links["F3H"][0],
            "DFR_accession_KT_registered":dfr[0] if dfr else None,
            "DFR_KT_registration_status":"RECORDED" if dfr else "NO_DFR_IN_KT_RANGE_NOT_GENE_ABSENCE",
            "DFR_CDS_annotation_in_KT_record":(None if not dfr else dfr[0] not in missingCDS),
            "functional_joint_toolkit":"UNMEASURED",
            "direct_flower_reexpression":"UNMEASURED",
            "historical_regain_evidence":("PUBLISHED_PHYLOGENETIC_INFERENCE"
                if taxon=="Iochroma edule" else "NOT_REESTIMATED")
        })
    counts={}
    for status_value in ("UNPIGMENTED","PIGMENTED"):
        group=[x for x in table if x["trait_state_original_fig1"]==status_value]
        counts[status_value]={
            "species":len(group),
            "CHI_KT_records":sum(bool(x["CHI_accession"]) for x in group),
            "F3H_KT_records":sum(bool(x["F3H_accession"]) for x in group),
            "DFR_KT_records":sum(bool(x["DFR_accession_KT_registered"]) for x in group),
            "DFR_KT_records_with_annotated_CDS":sum(
                x["DFR_CDS_annotation_in_KT_record"] is True for x in group),
            "DFR_KT_records_without_CDS_feature":sum(
                x["DFR_CDS_annotation_in_KT_record"] is False for x in group),
            "DFR_no_KT_deposit_species":sum(
                x["DFR_accession_KT_registered"] is None for x in group)
        }
    if (counts["UNPIGMENTED"]["DFR_KT_records"]!=5
        or counts["PIGMENTED"]["DFR_KT_records"]!=9
        or counts["UNPIGMENTED"]["DFR_KT_records_without_CDS_feature"]!=2
        or counts["PIGMENTED"]["DFR_KT_records_without_CDS_feature"]!=0):
        raise ValueError("Original accession and trait joined denominators changed")
    ie=next(x for x in table if x["current_taxon_name"]=="Iochroma edule")
    if not (ie["trait_state_original_fig1"]=="PIGMENTED"
            and ie["DFR_accession_KT_registered"]=="KT898445"):
        raise ValueError("Iochroma edule species-level inference bridge no longer source backed")
    return {
        "schema":"anto_iochrominae_2016_flower_anthocyanin_KT_coding_source_join_v1",
        "source_study_DOI":"10.1186/s12862-016-0675-3",
        "data_unit":"original named plant taxon × registered gene accession",
        "not_independent_evolutionary_transition_events":True,
        "n_species":len(table),
        "original_source_counts":counts,
        "known_source_coding_annotation_bias_possible":True,
        "DFR_missing_accession_not_gene_loss":True,
        "CDS_feature_missing_not_enzyme_loss":True,
        "published_regain_in_Iochroma_edule_not_independent_result":True,
        "original_paper_alignment_complete_but_not_reproduced":True,
        "functional_biochemistry_or_expression_verified_here":False,
        "independent_cross_clade_regain_rate_estimated":False,
        "taxon_gene_records":table
    }

def main():
    import argparse
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path,required=True)
    a.add_argument("--expected",type=Path)
    args=a.parse_args()
    result=calculate(load(GENES),load(TRAITS),load(ANNOTATIONS))
    if args.expected and load(args.expected)!=result:
        raise RuntimeError("Original-source frozen joined outcome drift")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("IOCHROMINAE_PIGMENT_SOURCE_JOIN_PASS",json.dumps(result["original_source_counts"],sort_keys=True))

if __name__=="__main__":main()
