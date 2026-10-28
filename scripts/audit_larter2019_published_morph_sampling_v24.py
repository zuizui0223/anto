"""Original-published Larter 2019 Table 1 ascertainment audit (ANTO v24).

This is independently recalculated, source-anchored *descriptive sampling
geometry* of already published Table 1 values, not a replication of the
Dryad-archived raw individual HPLC/qPCR data. The registered Dryad original
RAR is still HTTP 401/403 blocked (ANTO v23).

Original 2019 table entries include PURPLE/YELLOW for I. parvifolium.
Both white and yellow are called "acyanic" here. The percentages are
published approximate ESTIMATES from combined herbarium and natural field
surveys, NOT random unbiased ecological frequencies or original white-count
numerators. Do not treat surveyed counts as binomial trials.

2019's intentionally enriched acyanic bud samples are for investigating
molecular/developmental contrasts, not for estimating prevalence. Never use
their 2019 recruitment composition as natural morph frequencies, or fit
multigeneration evolutionary accessibility from these six taxa.

Five taxa also occur in 2016 Ho-Smith original-source 22-taxon list, but
these are NOT the same vouchered plants, assays or independent loss events.
Saracha punctata is NOT Saracha quitensis or S. nigribaccata.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL_TABLE=ROOT/"data/external/larter_2019_published_table1_morph_sampling_v24.json"
SOURCE2016=ROOT/"results/iochrominae_2016_trait_KT_source_join_20261019.json"
SOURCE2019_GATE=ROOT/"results/iochroma_2019_source_access_v23_20261028.json"

FROZEN_ROWS=(
    ("Eriolarynx fasciculata","purple/white",3,3,"1.8",225),
    ("Iochroma calycinum","purple/white",3,1,"0.9",225),
    ("Iochroma cyaneum","purple/white",3,2,"1.3",150),
    ("Iochroma parvifolium","purple/yellow",3,1,"1.3",75),
    ("Iochroma umbellatum","purple/white",3,3,"57",200),
    ("Saracha punctata","purple/white",3,2,"33",150),
)


def read(path:Path):
    return json.loads(path.read_text(encoding="utf-8"))


def source_rows(table):
    """Fail closed on species identity, original observed and estimated cells."""
    if (table.get("schema")!="anto_larter_2019_original_published_table1_acyanic_morph_sampling_v24"
        or table.get("paper_DOI")!="10.1002/dvdy.82"
        or table.get("source_type")!="published_2019_article_Table_1_manual_source_transcription_not_Dryad_RAR"
        or table.get("source_PDF")!="https://par.nsf.gov/servlets/purl/10431252"
        or table.get("original_sample_pair_order")!="pigmented/acyanic"):
        raise ValueError("Source DOI and original published Table 1 provenance drift")
    if not (table.get("raw_2019_individual_HPLC_qPCR_retrieved") is False
            and table.get("original_2019_RAR_bytes_retrieved") is False
            and table.get("new_plant_assays_or_regain_events") is False
            and table.get("no_original_survey_acyanic_counts_inferred_from_rounded_percentages") is True
            and table.get("no_binomial_convenience_interval_from_mixed_survey") is True):
        raise ValueError("Do not promote published Table 1 into unseen original individual data")
    rows=table.get("rows")
    if not isinstance(rows,list) or len(rows)!=6:
        raise ValueError("Expected 6 published polymorphic species")
    out=[]
    for entry,expected in zip(rows,FROZEN_ROWS):
        obs=(entry.get("published_taxon"),entry.get("color_labels"),
             entry.get("pigmented_sampled"),entry.get("acyanic_sampled"),
             entry.get("estimated_acyanic_frequency_percent"),
             entry.get("surveyed_individuals"))
        if obs!=expected or any(type(x) is not int for x in
               (entry["pigmented_sampled"],entry["acyanic_sampled"],entry["surveyed_individuals"])):
            raise ValueError("Original 2019 Table 1 species, color or numeric cells changed")
        p=F(entry["estimated_acyanic_frequency_percent"])/100
        n=entry["pigmented_sampled"]+entry["acyanic_sampled"]
        q=F(entry["acyanic_sampled"],n)
        if not 0<p<1 or not 0<q<1 or n>=entry["surveyed_individuals"]:
            raise ValueError("Invalid original morphological sampling range")
        out.append((entry,p,q,n))
    if sum(n for _,_,_,n in out)!=30:
        raise ValueError("Original 30 deliberately enrolled morph specimens drift")
    return out


def rounding_delta(percent_string):
    """Half the lowest originally printed percentage unit.

    E.g. 1.8% -> ±0.05 percentage POINT; 57% -> ±0.5 point.
    These are arithmetic rounding sensitivities, NOT sampling CIs.
    """
    unit=F(1,10**len(percent_string.split(".")[1])) if "." in percent_string else F(1)
    return unit/F(2*100)


def group_summary(rows):
    n_sample=sum(n for _,_,_,n in rows)
    w_sample=sum(x["acyanic_sampled"] for x,_,_,_ in rows)
    n_survey=sum(x["surveyed_individuals"] for x,_,_,_ in rows)
    p_weighted=sum(x["surveyed_individuals"]*p for x,p,_,_ in rows)/n_survey
    sample_share=F(w_sample,n_sample)
    low=sum(x["surveyed_individuals"]*(p-rounding_delta(
        x["estimated_acyanic_frequency_percent"])) for x,p,_,_ in rows)/n_survey
    high=sum(x["surveyed_individuals"]*(p+rounding_delta(
        x["estimated_acyanic_frequency_percent"])) for x,p,_,_ in rows)/n_survey
    if low<=0:
        raise ValueError("Rounding sensitivity invalid")
    return {
        "number_of_published_species":len(rows),
        "intentionally_recruited_pigmented":n_sample-w_sample,
        "intentionally_recruited_acyanic":w_sample,
        "intentionally_recruited_total":n_sample,
        "original_mixed_field_herbarium_surveyed_total":n_survey,
        "experimental_acyanic_sample_fraction":sample_share,
        "survey_n_weighted_original_PUBLISHED_estimated_frequency":p_weighted,
        "relative_acyanic_enrichment_of_experimental_composition":sample_share/p_weighted,
        "relative_enrichment_rounding_interval_not_sampling_CI":[sample_share/high,sample_share/low],
        "not_a_census_or_sampling_based_population_prevalence":True,
    }


def audit(original2019, source2016, source_gate):
    pairs=source_rows(original2019)
    if (source_gate.get("original_public_dataset_DOI")!="10.5061/dryad.p5dq84v"
        or source_gate.get("original_RAR_bytes_verified") is not False
        or source_gate.get("raw_gene_expression_sample_rows_parsed")!=0):
        raise ValueError("Original 2019 Dryad source gate changed; raw reanalysis requires new review")
    if source2016.get("n_species")!=22 or len(source2016.get("taxon_gene_records",[]))!=22:
        raise ValueError("2016 original GenBank crosswalk must contain 22 source taxa")
    known={x["current_taxon_name"]:x for x in source2016["taxon_gene_records"]}
    if len(known)!=22:
        raise ValueError("2016 crosswalk has duplicate source taxa")
    results=[]
    for x,p,q,n in pairs:
        name=x["published_taxon"]
        cross=known.get(name)
        if cross is not None and (cross["trait_state_original_fig1"]!="PIGMENTED"
                   or cross["functional_joint_toolkit"]!="UNMEASURED"
                   or cross["direct_flower_reexpression"]!="UNMEASURED"):
            raise ValueError("2016 recorded color or competence is not source-matched")
        delta=rounding_delta(x["estimated_acyanic_frequency_percent"])
        results.append({
            "published_2019_taxon":name,
            "acyanic_morph_is_yellow_not_white":x["color_labels"]=="purple/yellow",
            "original_2019_morph_labels":x["color_labels"],
            "original_morph_recruited_pigmented":x["pigmented_sampled"],
            "original_morph_recruited_acyanic":x["acyanic_sampled"],
            "experimental_morph_total":n,
            "original_estimated_acyanic_frequency_percent":x["estimated_acyanic_frequency_percent"],
            "original_herbarium_and_field_individuals_surveyed":x["surveyed_individuals"],
            "published_morph_frequency_as_fraction":p,
            "experimental_acyanic_fraction":q,
            "experimental_fraction_to_survey_frequency_ratio":q/p,
            "rounding_relative_enrichment_min":q/(p+delta),
            "rounding_relative_enrichment_max":q/(p-delta),
            "prob_random_at_least_one_acyanic_under_ARTIFICIAL_IID_sampling":
                1-(1-p)**n,
            "overlap_original_2016_GenBank_taxon_name":cross is not None,
            "2016_original_DFR_new_KT_accession_if_same_species":(
                cross["DFR_accession_KT_registered"] if cross else None),
            "2016_registrations_not_same_plant_or_2019_flower_function":True,
        })
    if results[-1]["overlap_original_2016_GenBank_taxon_name"]:
        raise ValueError("Do not confuse Saracha punctata with other Saracha taxa")
    if sum(r["overlap_original_2016_GenBank_taxon_name"] for r in results)!=5:
        raise ValueError("Cross-study species overlap must be 5 of 6, not 6/6")
    rare=[t for t in pairs if t[1]<F(2,100)]
    common=[t for t in pairs if t[1]>=F(30,100)]
    if len(rare)!=4 or len(common)!=2 or len(rare)+len(common)!=6:
        raise ValueError("Original Table 1 4-rare/2-common pattern drift")
    rare_stats=group_summary(rare)
    all_stats=group_summary(pairs)
    if rare_stats["intentionally_recruited_acyanic"]!=7 or rare_stats[
            "original_mixed_field_herbarium_surveyed_total"]!=675 or rare_stats[
            "survey_n_weighted_original_PUBLISHED_estimated_frequency"]!=F(1,75):
        raise ArithmeticError("Published 2019 rare morph descriptive denominators changed")
    if rare_stats["relative_acyanic_enrichment_of_experimental_composition"]!=F(525,19):
        raise ArithmeticError("Source-derived sampling enrichment exact fraction drift")
    if all_stats["intentionally_recruited_acyanic"]!=12 or all_stats[
            "original_mixed_field_herbarium_surveyed_total"]!=1025:
        raise ArithmeticError("Published 2019 total recruitment and mixed survey counts drift")
    if [r["published_2019_taxon"] for r in results if r[
            "original_morph_recruited_acyanic"]==1]!=[
            "Iochroma calycinum","Iochroma parvifolium"]:
        raise ArithmeticError("Critical one-acyanic-morph biological replicate gate changed")
    return {
        "schema":"anto_v24_larter2019_published_table1_acyanic_morph_ascertainment",
        "date":"2026-10-28",
        "original_2019_article_DOI":"10.1002/dvdy.82",
        "original_2019_table":"Table 1, original published page 1094",
        "underlying_source_method":"Human checked transcription of PUBLIC 2019 Table 1, not inaccessible original Dryad RAR",
        "published_source_table":"data/external/larter_2019_published_table1_morph_sampling_v24.json",
        "source_full_text_PDF":"https://par.nsf.gov/servlets/purl/10431252",
        "individually_measured_2019_gene_expression_or_HPLC_obtained":False,
        "actual_2019_acyanic_morph_enrichment_is_purposive_recruitment_not_a_population_frequency":True,
        "original_survey_mixed_herbarium_and_field_not_unbiased_iid":True,
        "source_percentages_are_rounded_estimates_not_integer_raw_white_counts":True,
        "single_acyanic_morph_recruited_in_two_original_species":True,
        "rar_original_source_gate_currently_blocked":True,
        "purple_yellow_acyanic_counted_but_not_literally_white":True,
        "published_species":results,
        "rare_below_2percent_published_frequency_4_species":rare_stats,
        "common_at_least_30percent_published_frequency_2_species":group_summary(common),
        "all_six_published_polymorphic_species":all_stats,
        "morph_presence_in_four_to_six_iid_random_samples_is_counterfactual_only":True,
        "2016_original_22taxon_new_KT_panel_overlap_species":5,
        "2016_unmatched_source_taxon":"Saracha punctata",
        "not_independent_2016_genomic_plus_2019_flower_assay_same_individual":True,
        "no_new_estimated_genetic_ancestral_regain_or_functional_retention_rate":True,
        "no_actual_2019_expression_regression_recomputed":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
    }


def serialize(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {k:serialize(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [serialize(v) for v in value]
    return value


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    parser.add_argument("--expected",type=Path)
    args=parser.parse_args()
    value=serialize(audit(read(ORIGINAL_TABLE),read(SOURCE2016),read(SOURCE2019_GATE)))
    if args.expected and value!=read(args.expected):
        raise RuntimeError("v24 original published 2019 Table 1 source arithmetic drift")
    result=json.dumps(value,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(result,encoding="utf-8")
    else:
        print(result)


if __name__=="__main__":
    main()
