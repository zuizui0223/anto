"""ANTO v19: Original-source evidence admissibility, not trait-regain fitting.

Use frozen source-verified records from *distinct observation designs*:
- Ho & Smith 2016: 22 original Iochrominae taxa, 58 original KT source
  accessions with phenotype and DFR source/annotation status.
- Wheeler et al. 2025: 24 blue/red Iochroma backcross floral RNA libraries.
- Wessinger & Rausher 2014: 161 Penstemon hybrid F2 phenotypes.
- v11 mechanism registry: selected Petunia, Iochrominae, Antirrhinum,
  Mimulus and Cirsium mechanistic evidence cases.

The endpoint is whether these original source masks justify fitting an
intrinsic floral-trait regain and long-term functional-capacity rate.
They do not. This is a reproducible ORIGINAL-SOURCE-RECORD audit, not
our replication of 2016 phylogenetic codon-selection or floral enzyme
experiments, not a species-level meta-analysis, and not a gene-loss claim.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction as F
from math import comb
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE_JOIN=ROOT/"results/iochrominae_2016_trait_KT_source_join_20261019.json"
REGISTRY=ROOT/"data/afterlife_v11_mechanistic_regain_routes_literature_registry_20261020.json"
PENSTEMON=ROOT/"data/external/penstemon_2014_f3p5h_source_genetic_accessibility_v1.json"
IOCHROMA=ROOT/"results/iochroma_2025_color_conditioned_sensitivity_v3_20261026.json"
EXPECTED_ORIGINAL_DOI="10.1186/s12862-016-0675-3"
STATUS=("DFR_KT_CDS_ANNOTATED","DFR_KT_NO_CDS_FEATURE","DFR_NOT_IN_NEW_KT_BLOCK")
COLOURS=("UNPIGMENTED","PIGMENTED")


def read(path:Path):
    return json.loads(path.read_text(encoding="utf-8"))


def record_status(row):
    """Distinguish never deposited in this new KT block from feature absent."""
    a=row.get("DFR_accession_KT_registered")
    b=row.get("DFR_CDS_annotation_in_KT_record")
    registration=row.get("DFR_KT_registration_status")
    if a is None:
        if b is not None or registration!="NO_DFR_IN_KT_RANGE_NOT_GENE_ABSENCE":
            raise ValueError("Missing DFR accession cannot be a gene deletion or CDS call")
        return STATUS[2]
    if not isinstance(a,str) or not a.startswith("KT8984"):
        raise ValueError("Wrong deposited DFR accession")
    if registration!="RECORDED":
        raise ValueError("Registered DFR accession must retain its original status")
    if b is True:
        return STATUS[0]
    if b is False:
        return STATUS[1]
    raise ValueError("Present DFR accession with unknown CDS metadata is inconsistent")


def contingency(join):
    if join.get("schema")!="anto_iochrominae_2016_flower_anthocyanin_KT_coding_source_join_v1":
        raise ValueError("Original phenotype × accession source schema mismatch")
    if join.get("source_study_DOI")!=EXPECTED_ORIGINAL_DOI:
        raise ValueError("Original Ho–Smith source DOI mismatch")
    rows=join.get("taxon_gene_records")
    if not isinstance(rows,list) or len(rows)!=22:
        raise ValueError("Exactly 22 original Iochrominae taxa required")
    names=[r["current_taxon_name"] for r in rows]
    if len(set(names))!=22:
        raise ValueError("Original taxon names must be distinct")
    types=Counter(r["taxonomic_join_type"] for r in rows)
    if types!={"DIRECT_BINOMIAL":17,"POST2016_COMBINATION_OR_SYNONYM":5}:
        raise ValueError("Original taxon synonyms or author name coverage drift")
    table={c:{s:0 for s in STATUS} for c in COLOURS}
    gene_accessions=set()
    without_cds=[]
    for r in rows:
        phenotype=r["trait_state_original_fig1"]
        if phenotype not in COLOURS:
            raise ValueError("Unexpected original flower pigment category")
        if r.get("functional_joint_toolkit")!="UNMEASURED" or r.get("direct_flower_reexpression")!="UNMEASURED":
            raise ValueError("Historical source record is not validated functional reexpression")
        for gene in ("CHI_accession","F3H_accession"):
            accession=r[gene]
            if not isinstance(accession,str) or not accession.startswith("KT898"):
                raise ValueError("Original CHI/F3H source accession missing")
            if accession in gene_accessions:
                raise ValueError("Source accession must be unique, not a replicate")
            gene_accessions.add(accession)
        if r["DFR_accession_KT_registered"] is not None:
            if r["DFR_accession_KT_registered"] in gene_accessions:
                raise ValueError("DFR original accession duplicated")
            gene_accessions.add(r["DFR_accession_KT_registered"])
        key=record_status(r)
        table[phenotype][key]+=1
        if key==STATUS[1]:
            without_cds.append((r["current_taxon_name"],r["DFR_accession_KT_registered"]))
    if len(gene_accessions)!=58:
        raise ValueError("Original 58 accession biological sequence coverage drift")
    if set(gene_accessions)!={f"KT{i}" for i in range(898394,898452)}:
        raise ValueError("Original KT898394–KT898451 source identities drift")
    if table!={"UNPIGMENTED":{
        STATUS[0]:3,STATUS[1]:2,STATUS[2]:2},
        "PIGMENTED":{
        STATUS[0]:9,STATUS[1]:0,STATUS[2]:6}}:
        raise ValueError("Original source 7/15 flower and 3-state DFR mask changed")
    if sorted(without_cds)!=[
            ("Iochroma confertiflorum","KT898443"),
            ("Iochroma ellipticum","KT898446")]:
        raise ValueError("The two missing CDS features must preserve source identities")
    if any(sum(table[group].values())!=target for group,target in
           (("UNPIGMENTED",7),("PIGMENTED",15))):
        raise ValueError("Original published taxon counts changed")
    if not (join["DFR_missing_accession_not_gene_loss"]
            and join["CDS_feature_missing_not_enzyme_loss"]
            and join["not_independent_evolutionary_transition_events"]
            and not join["functional_biochemistry_or_expression_verified_here"]
            and not join["independent_cross_clade_regain_rate_estimated"]):
        raise ValueError("Functional, accession and macroevolution claim boundaries changed")
    return table


def x2(table):
    """Exact Pearson X2 of a *source metadata* 2×3 table, not genetic loss."""
    first=[table["UNPIGMENTED"][s] for s in STATUS]
    cols=[sum(table[p][s] for p in COLOURS) for s in STATUS]
    n=sum(cols)
    n_white=sum(first)
    n_pig=n-n_white
    return sum(
        (F(white)-F(n_white*col,n))**2/F(n_white*col,n)
        +(F(col-white)-F(n_pig*col,n))**2/F(n_pig*col,n)
        for white,col in zip(first,cols)
    )


def conditional_metadata_reference(table):
    """Enumerate all 2×3 tables at fixed margins, exact hypergeometric mass.

    Exchangeable taxon labels are *not biologically defensible* here:
    lineages share phylogenetic ancestry and author sampling is not a
    random experiment. Thus this is an arithmetic reference ONLY, and the
    p value must never be interpreted as macroevolutionary significance.
    """
    row=tuple(table["UNPIGMENTED"][s] for s in STATUS)
    cols=[sum(table[p][s] for p in COLOURS) for s in STATUS]
    n=sum(cols)
    n_white=sum(row)
    obs=x2(table)
    tail=F(0);mass=F(0);num=0
    for a in range(cols[0]+1):
        for b in range(cols[1]+1):
            c=n_white-a-b
            if not 0<=c<=cols[2]:
                continue
            num+=1
            prob=F(comb(cols[0],a)*comb(cols[1],b)*comb(cols[2],c),
                   comb(n,n_white))
            mass+=prob
            tab={COLOURS[0]:dict(zip(STATUS,(a,b,c))),
                 COLOURS[1]:dict(zip(STATUS,
                    (cols[0]-a,cols[1]-b,cols[2]-c)))}
            if x2(tab)>=obs:
                tail+=prob
    if mass!=1 or num!=21:
        raise ArithmeticError("Exact metadata table enumeration failed")
    return {"metadata_observed_pearson_X2_exact":str(obs),
            "fixed_margins_conditional_tail_exact":str(tail),
            "enumerated_2_by_3_tables":num,
            "conditional_p_is_not_a_phylogenetic_or_evolutionary_test":True}


def analysis(join, literature, penstemon, iochroma):
    table=contingency(join)
    un=table["UNPIGMENTED"]
    pi=table["PIGMENTED"]
    source_verified={
        "original_taxa":22,
        "source_nucleotide_accessions":58,
        "unpigmented_taxa":7,
        "pigmented_taxa":15,
        "three_source_metadata_levels":{k:v for k,v in table.items()},
        "CHI_F3H_new_KT_accessions_by_taxon":22,
        "DFR_new_KT_accessions":14,
        "DFR_CDS_features_annotated_in_new_block":12,
        "DFR_records_with_no_CDS_feature":2,
        "DFR_no_new_KT_accession_species":8,
        "new_KT_DFR_records_with_CDS_feature_unpigmented":3,
        "new_KT_DFR_records_with_CDS_feature_pigmented":9,
        "observed_source_annotation_fraction_unpigmented":str(F(3,7)),
        "observed_source_annotation_fraction_pigmented":str(F(9,15)),
        "observed_annotation_fraction_difference_unpig_minus_pig":str(F(3,7)-F(9,15)),
        "unpigmented_without_new_KT_annotated_DFR":un[STATUS[1]]+un[STATUS[2]],
        "unpigmented_missing_KT_annotated_DFR_fraction":str(F(4,7)),
        "original_author_reports_intact_three_coding_regions_after_floral_pigment_loss":True,
        "author_full_coding_in_original_alignments_not_independently_recomputed":True,
        "false_negative_if_missing_new_DFR_CDS_feature_treated_as_gene_absence_unpigmented":4,
        "unknown_DFR_new_KT_feature_statuses_do_not_establish_pseudogenization":True,
    }
    # These are *hypothetical* completions of a biased source-annotation
    # *proxy*, only after explicitly assuming an annotated record counts as
    # positive evidence and missing records are of unknown status. Neither
    # premise provides direct flowering-pathway functional measurements.
    bound_low=F(3,7)-F(15,15)
    bound_high=F(7,7)-F(9,15)
    if bound_low>=0 or bound_high<=0:
        raise ArithmeticError("Must include both signs under arbitrary proxy completion")
    partial={
        "counterfactual_proxy_assumption":"All annotated DFR CDS new source records provisionally counted 1; unverified rows can be 0 or 1, not actual gene-function calls",
        "unpigmented_proxy_completion_fraction_range":["3/7","1"],
        "pigmented_proxy_completion_fraction_range":["3/5","1"],
        "unpig_minus_pig_proxy_contrast_sharp_bounds":[str(bound_low),str(bound_high)],
        "source_proxy_contrast_sign_unidentified":True,
        "true_functional_joint_toolkit_unmeasured_in_all_22_by_THIS_KT_SOURCE_AUDIT":True,
        "true_functional_joint_toolkit_fraction_identification_bounds_both_groups":["0","1"],
        "true_group_function_difference_undetermined_without_functional_assay":True,
        "no_gene_deletion_or_pseudogene_inferred":True,
    }
    dfr_annotated_un=un[STATUS[0]]
    dfr_no_annotation_un=un[STATUS[1]]
    dfr_annotated_pig=pi[STATUS[0]]
    dfr_no_annotation_pig=pi[STATUS[1]]
    if [dfr_annotated_un,dfr_no_annotation_un,dfr_annotated_pig,dfr_no_annotation_pig]!=[3,2,9,0]:
        raise ValueError("Recorded-Dfr annotation conditional denominator changed")
    # One-sided hypergeometric probability of both 2 unannotated DFR
    # source records falling among 5 unpigmented-record taxa at fixed
    # 5/9 accession and 2 unannotated margins.
    one_sided=F(comb(5,2),comb(14,2))
    reference=conditional_metadata_reference(table)
    reference["among_14_deposited_DFR_two_no_CDS_both_unpigmented_exact"]=str(one_sided)

    cases=literature.get("cases")
    if not isinstance(cases,list) or len(cases)!=6:
        raise ValueError("v11 source-case list changed; requires explicit review")
    ids=[r.get("case_id") for r in cases]
    if len(set(ids))!=6:
        raise ValueError("v11 must use distinct original mechanism-case IDs")
    petunia=[x for x in cases if x.get("clade")=="Petunia"]
    if len(petunia)!=2 or sum(x.get("molecular_causal_evidence","").startswith("SOURCE_SUPPORTED") for x in petunia)!=2:
        raise ValueError("Two source-supported Petunia cases from ONE genus required")
    if any(x.get("independent_regain_clade_weight") not in (0,1) for x in cases):
        raise ValueError("Invalid clade weights")
    if sum(x["independent_regain_clade_weight"] for x in cases)!=1:
        raise ValueError("Petunia is one independent molecular case clade, not two")
    if penstemon.get("published_F2_n")!=161:
        raise ValueError("Original Penstemon F2 experimental cross 161 count changed")
    counts=penstemon["observed_F2_table"]
    if (sum(sum(row.values()) for row in counts.values())!=161 or
        counts["purple"]!={"FF":57,"Ff":65,"ff":0} or
        counts["red"]!={"FF":0,"Ff":15,"ff":24}):
        raise ValueError("Original Penstemon genotype × phenotype matrix drift")
    if not penstemon["reported_source_molecular_results"]["no_direct_natural_pbarbatus_backmutation_reversal_observed"]:
        raise ValueError("Penstemon donor complementation must not be natural backmutation")
    if iochroma["original_24_sensitivity"]["original_individuals"]!=24:
        raise ValueError("Original Iochroma backcross 24-unit source altered")
    if iochroma["original_24_sensitivity"]["other_preselected_transcript_entries"]!=15:
        raise ValueError("Original Iochroma transcript comparator count drift")
    if not iochroma["not_independent_plant_clade_or_natural_flower_regain"]:
        raise ValueError("Iochroma red/blue association is not a natural flower regain")
    if literature["derived_summary_contract"]["source_cases_with_trait_specific_qualified_viable_seed_fitness"]!=0:
        raise ValueError("Source registry has changed fitness evidence status")

    cross={
        "source_programmes":[
            {"programme":"HoSmith2016_Iochrominae","observations":22,
             "observational_unit":"species with colour symbol and deposited-source accessions",
             "endpoint":"source data presence/annotation, not coding enzyme function or independent regain event",
             "clade":"Iochrominae"},
            {"programme":"Wheeler2025_Iochroma","observations":24,
             "observational_unit":"flower RNA sample from one interspecies backcross offspring",
             "endpoint":"in-sample red/blue transcript differences, not historical white-flower regain",
             "clade":"Iochrominae"},
            {"programme":"WessingerRausher2014_Penstemon","observations":161,
             "observational_unit":"experimental F2 hybrid offspring",
             "endpoint":"blue-purple donor allele complementation, not natural intrinsic backmutation",
             "clade":"Penstemon"},
        ],
        "incompatible_original_observation_units_raw_arithmetic_sum":207,
        "raw_sum_must_not_be_called_regain_events_or_independent_replicates":True,
        "petunia_molecular_mechanism_source_cases":2,
        "petunia_independent_genus_units":1,
        "independent_evolutionary_regain_clades_molecularly_source_supported_in_selected_v11_registry":1,
        "selected_registry_inferred_iochroma_regain_is_not_a_causally_proven_second_clade":True,
        "multiple_iochroma_studies_are_within_same_iochrominae_lineage_group":True,
        "source_cases_with_trait_specific_qualified_viable_seed_fitness":0,
        "observed_species_level_phylogenetic_transition_directly_reconstructed_by_anto":False,
        "real_plant_generational_ecological_rate_estimated_by_anto":False,
    }
    return {
        "schema":"anto_v19_original_sources_cross_system_observability_gate",
        "date":"2026-10-26",
        "source_reuse_not_original_author_experiment_or_original_full_alignment_reconstruction":True,
        "source_paths":{
            "phenotype_genbank_22_taxa":"results/iochrominae_2016_trait_KT_source_join_20261019.json",
            "selected_mechanisms":"data/afterlife_v11_mechanistic_regain_routes_literature_registry_20261020.json",
            "penstemon_hybrid_F2":"data/external/penstemon_2014_f3p5h_source_genetic_accessibility_v1.json",
            "iochroma_flower_rna":"results/iochroma_2025_color_conditioned_sensitivity_v3_20261026.json"
        },
        "Iochrominae_KT_source_metadata":source_verified,
        "exact_source_metadata_random_label_reference_NOT_evolutionary_p_value":reference,
        "unmeasured_function_partial_identification":partial,
        "cross_programme_evidence_grain":cross,
        "no_natural_floral_trait_regain_rate_estimated":True,
        "no_long_term_fitness_or_adaptation_identified":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    parser.add_argument("--expected",type=Path)
    args=parser.parse_args()
    result=analysis(read(SOURCE_JOIN),read(REGISTRY),read(PENSTEMON),read(IOCHROMA))
    if args.expected and result!=read(args.expected):
        raise RuntimeError("Original-source v19 frozen evidence admissibility drift")
    value=json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(value,encoding="utf-8")
    else:
        print(value)


if __name__=="__main__":
    main()
