"""Source-claim validation of empirical (non-random) trait-regain cases.

This registry evaluates provenance, comparative units and which published
regain outcomes are *eligible for description*. It cannot estimate frequencies,
rate parameters or historical fitness from a biased small case selection.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/afterlife_v11_mechanistic_regain_routes_literature_registry_20261020.json"

def validate(d):
    assert d.get("schema")=="anto_evolutionary_afterlife_empirical_route_registry_v11"
    cases=d["cases"]
    if len(cases)!=6 or len({c["case_id"] for c in cases})!=6:
        raise ValueError("Exactly six uniquely documented case records required")
    if any(not c.get("source_dois") or c["trait_specific_viable_seed_fitness"]=="DIRECT_MEASURED"
           for c in cases):
        raise ValueError("Missing source or unwarranted ecological fitness claim")
    classes=Counter(c["clade"] for c in cases)
    if len(classes)!=5 or classes["Petunia"]!=2 or any(n!=1 for cl,n in classes.items() if cl!="Petunia"):
        raise ValueError("Clade provenance is not an independent genus tally")
    sourced=[c for c in cases if c["case_role"] in (
        "MOLECULARLY_SUPPORTED_PIGMENT_REGAIN",
        "MOLECULARLY_SUPPORTED_ALTERNATIVE_HUE_REGAIN",
        "PHYLOGENETIC_PIGMENT_REGAIN_INFERRED")]
    if len(sourced)!=3:
        raise ValueError("Do not recode loss and pattern regulation as floral pigment evolutionary regains")
    restored=[c for c in sourced if c["ancestral_hue_restoration"]=="SOURCE_SUPPORTED"]
    altered=[c for c in sourced if c["ancestral_hue_restoration"]=="ALTERED_NOT_ANCESTRAL"]
    if len(restored)!=1 or len(altered)!=2:
        raise ValueError("Original hue vs altered hue source readout changed")
    if {(c["taxon"],c["clade"]) for c in restored}!={("Petunia secreta","Petunia")}:
        raise ValueError("Ancestral-hue restoration case changed")
    if {c["taxon"] for c in altered}!={"Petunia exserta","Iochroma edule"}:
        raise ValueError("Altered hue regain candidates drifted")
    causal=[c for c in sourced if c["molecular_causal_evidence"].startswith("SOURCE_SUPPORTED")]
    if {c["clade"] for c in causal}!={"Petunia"} or len(causal)!=2:
        raise ValueError("The single causal-recovery clade cannot be doubled")
    i=next(c for c in cases if c["taxon"]=="Iochroma edule")
    if i["molecular_causal_evidence"]!="NOT_IDENTIFIED" or i["ancestral_hue_restoration"]!="ALTERED_NOT_ANCESTRAL":
        raise ValueError("Iochroma causal ancestral-route status inaccurately promoted")
    ant=next(c for c in cases if c["clade"]=="Antirrhinum")
    if ant["case_role"]!="EXPERIMENTALLY_INHERITED_FLORAL_PATTERN_REGULATOR_REVERSIBILITY_CONTROL":
        raise ValueError("Antirrhinum pattern repressor mistaken as ancient floral regain")
    cir=next(c for c in cases if c["clade"]=="Cirsium")
    if (cir["molecular_causal_evidence"]!="NOT_MEASURED" or
        cir["ancestral_hue_restoration"]!="NOT_IDENTIFIABLE"):
        raise ValueError("Leaf-only Cirsium source mislabeled as restored evolutionary phenotype")
    if any(c["quantitative_lineage_reexpression_rate"]!="NOT_ESTIMATED" for c in cases):
        raise ValueError("Source case selection cannot estimate regain probability")
    report={
        "source_registry":"data/afterlife_v11_mechanistic_regain_routes_literature_registry_20261020.json",
        "status":"CURATED_PUBLISHED_CASE_COMPARISON_NOT_META_ANALYSIS",
        "n_cases":len(cases),
        "n_distinct_clades":len(classes),
        "regain_cases_with_published_ancestry":len(sourced),
        "ancestral_hue_restoration_cases":len(restored),
        "altered_hue_regain_cases":len(altered),
        "causal_regain_examples_same_Petunia_genus":len(causal),
        "causal_regain_independent_genus_count":len({c["clade"] for c in causal}),
        "phylogenetically_inferred_regain_in_separate_Iochrominae_clade":1,
        "measured_Cirsium_floral_reexpression_events":0,
        "legitimate_cross_clade_frequency_inference":False,
        "source_cases_are_not_6_independent_regains":True,
        "causal_reproductive_fitness_after_regain_verified_by_current_registry":False,
        "phenotypic_regain_cannot_be_equated_with_ancestral_hue_restoration":True,
        "registry_uses_cited_original_authors_results_not_new_raw_plant_experiment":True,
        "per_clade_case_counts":dict(sorted(classes.items())),
        "source_case_ids":sorted(c["case_id"] for c in cases)
    }
    expected=d["derived_summary_contract"]
    for key,k in [
        ("distinct_named_clades","n_distinct_clades"),
        ("listed_pigment_regain_candidates_with_sourced_evolutionary_ancestry","regain_cases_with_published_ancestry"),
        ("published_molecularly_supported_pigment_regains_within_petunia","causal_regain_examples_same_Petunia_genus"),
        ("independently_molecularly_proven_pigment_regain_clades_in_this_registry","causal_regain_independent_genus_count"),
        ("source_cases_currently_showing_ancestral_hue_restoration","ancestral_hue_restoration_cases"),
        ("source_cases_currently_showing_altered_hue_on_colour_return","altered_hue_regain_cases"),
    ]:
        if expected[key]!=report[k]:
            raise ValueError("Registry derived totals not consistent with case-level evidence")
    if expected["human_selected_source_cases_not_unbiased_clade_sample"] is not True:
        raise ValueError("Non-systematic inclusion gate was lost")
    return report

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    d=json.loads(SOURCE.read_text(encoding="utf-8"))
    report=validate(d)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("AFTERLIFE_V11_SOURCE_EVIDENCE_GATE_PASS",json.dumps(
        {k:v for k,v in report.items() if k in (
            "n_cases","n_distinct_clades","ancestral_hue_restoration_cases",
            "altered_hue_regain_cases","causal_regain_independent_genus_count")},sort_keys=True))

if __name__=="__main__":
    main()
