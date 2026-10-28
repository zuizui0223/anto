"""ANTO v25: source-published pathway targets, NOT new raw expression fits.

Original Larter et al. 2019 DOI 10.1002/dvdy.82 Table 2 lists only
the five best published OU-PGLS AIC candidates for two CONTINUOUS
pigment outputs. Original Berardi et al. 2016 DOI
10.1016/j.phytochem.2016.05.007 Table 2 lists a 9x9 phylogenetic
correlation matrix based on 26 individuals from nine species. These
papers share Iochrominae and are NOT independent regain clades, same
individuals, or joint long-term evolution trials.

Our new calculation consists ONLY of deterministic AIC evidence ratios
and conditional, five-CANDIDATE Akaike weights, exact source table
consistency, and descriptive tissue/metabolite correlations.

AIC weights are NOT gene-causal probabilities, posterior effects,
molecular pathway reactivation rates or held-out prediction scores.
No original 2019 HPLC/qPCR RAR bytes accessed: v23 official Dryad
source gate HTTP 401/403. Original authors already reported CHI as
best quantitative flux predictor and DFR/ANS as fixed-white species
expression differences. Do not count that biological conclusion as
a new ANTO discovery.
"""
from __future__ import annotations
import argparse
import json
import math
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/external/iochrominae_2016_tissues_2019_pigment_aic_published_v25.json"
RAW_GATE=ROOT/"results/iochroma_2019_source_access_v23_20261028.json"
V24_TABLE=ROOT/"results/larter2019_morph_ascertainment_v24_20261028.json"

EXPECTED_MODELS={
    "anthocyanidins":[("Chi","66.2","0.0",6),
                      ("F3h","78.4","12.2",6),
                      ("Chi+F3h","92.0","25.8",10),
                      ("Chs","94.4","28.2",6),
                      ("Ans","119.3","53.1",6)],
    "flavonols":[("Chi","65.0","0.0",6),
                 ("F3h","74.9","9.9",6),
                 ("Chi+F3h","90.6","25.6",10),
                 ("Chs","94.7","29.7",6),
                 ("Ans","128.0","63.0",6)],
}
SOURCE_GENES={
    "Chi":("CHI",),"F3h":("F3H",),"Chi+F3h":("CHI","F3H"),
    "Chs":("CHS",),"Ans":("ANS",)
}
EXPECTED_CORRELATIONS=(
    ("1","-0.15","0.16","-0.13","-0.18","0.25","-0.14","-0.04","-0.04"),
    ("-0.15","1","0.42","0.40","0.78","0.29","0.13","0.00","-0.19"),
    ("0.16","0.42","1","0.54","0.44","0.96","0.04","-0.05","-0.05"),
    ("-0.13","0.40","0.54","1","0.68","0.60","0.22","-0.08","-0.07"),
    ("-0.18","0.78","0.44","0.68","1","0.46","0.22","0.13","-0.20"),
    ("0.25","0.29","0.96","0.60","0.46","1","0.06","-0.03","-0.04"),
    ("-0.14","0.13","0.04","0.22","0.22","0.06","1","0.52","-0.17"),
    ("-0.04","0.00","-0.05","-0.08","0.13","-0.03","0.52","1","-0.15"),
    ("-0.04","-0.19","-0.05","-0.07","-0.20","-0.04","-0.17","-0.15","1"),
)
EXPECTED_DIMENSIONS=(
    "petal_mono_anthocyanidin","petal_mono_flavonol","leaf_mono_flavonol",
    "petal_di_anthocyanidin","petal_di_flavonol","leaf_di_flavonol",
    "petal_tri_anthocyanidin","petal_tri_flavonol","leaf_tri_anthocyanidin",
)


def load(path:Path):
    return json.loads(path.read_text(encoding="utf-8"))


def source_verification(s, source_gate, v24):
    """Strong paper and exact printed-table identity, without raw-data promotion."""
    if not (s.get("schema")=="anto_v25_original_published_2016_tissue_2019_gene_aic_endpoint_registry"
            and s.get("original_published_sources_not_new_assays") is True):
        raise ValueError("Original published source identity changed")
    late=s.get("source_2019",{})
    early=s.get("source_2016",{})
    boundary=s.get("scope_control",{})
    if not (late.get("doi")=="10.1002/dvdy.82"
            and late.get("original_table")=="Table 2, p. 1096"
            and late.get("relative_model_selection_only_not_predictive_holdout") is True
            and late.get("registered_original_RAR_access_still_unverified") is True
            and early.get("doi")=="10.1016/j.phytochem.2016.05.007"
            and early.get("original_individuals_in_published_phylogenetic_correlation")==26
            and early.get("original_species_in_published_phylogenetic_correlation")==9
            and early.get("reported_phylogenetic_lambda")=="0.342"
            and early.get("original_matrix_rounded_to_two_decimals") is True
            and early.get("sample_covariance_or_new_phylogenetic_pvalue_not_computed") is True):
        raise ValueError("Original Table2 paper identity, response grain or meta-evidence drift")
    groups=late.get("reported_models")
    if not isinstance(groups,list) or len(groups)!=2:
        raise ValueError("Only the two original 2019 quantitative pigment response tables permitted")
    for group in groups:
        ref=EXPECTED_MODELS.get(group.get("response"))
        if ref is None or not isinstance(group.get("models"),list):
            raise ValueError("2019 pigment response not original Table2")
        seen=[]
        for row in group["models"]:
            if type(row.get("npar")) is not int:
                raise ValueError("Source AIC parameter count not an integer")
            name=row.get("model")
            genes=row.get("genes")
            if name not in SOURCE_GENES or genes!=list(SOURCE_GENES[name]):
                raise ValueError("AIC model gene identity drift")
            seen.append((name,row.get("aic"),row.get("delta_aic"),row.get("npar")))
        if seen!=ref:
            raise ValueError("2019 printed PGLS Table2 AIC or model label altered")
    if [g["response"] for g in groups]!=["anthocyanidins","flavonols"]:
        raise ValueError("2019 pigment response order/multiplicity altered")
    corr=early.get("original_printed_correlation_matrix")
    if tuple(tuple(x for x in row) for row in corr)!=EXPECTED_CORRELATIONS or (
        tuple(early.get("dimension_labels",()))!=EXPECTED_DIMENSIONS):
        raise ValueError("2016 original 9x9 tissue correlation identity changed")
    if any(not all(row[j]==corr[j][i] for j in range(9))
           for i,row in enumerate(corr)):
        raise ValueError("2016 published symmetric matrix altered")
    if any(corr[i][i]!="1" for i in range(9)):
        raise ValueError("2016 reported self-correlation diagonal altered")
    if not all(boundary.get(k) is True for k in (
            "same_clade_2016_and_2019_not_two_independent_reversions",
            "no_original_2019_hplc_qpcr_raw_reads",
            "no_new_direct_flower_joint_functional_assays",
            "no_natural_ancestral_colour_regain_rate",
            "no_2019_held_out_prediction_reestimated",
            "no_between_source_causal_transport")):
        raise ValueError("Biological data status must not be promoted")
    if boundary.get("Cirsium_empirical_owner")!="zuizui0223/aza3":
        raise ValueError("Cirsium organismal data custody mismatch")
    if not (source_gate.get("original_RAR_bytes_verified") is False
            and source_gate.get("raw_gene_expression_sample_rows_parsed")==0
            and source_gate.get("original_pub_DOI")=="10.1002/dvdy.82"
            and v24.get("original_2019_article_DOI")=="10.1002/dvdy.82"
            and v24.get("individually_measured_2019_gene_expression_or_HPLC_obtained") is False):
        raise ValueError("Original raw-material availability changed: new raw analysis required")
    return True


def conditional_akaike(group):
    """Relative weights strictly WITHIN five author-listed top model candidates.

    The full candidate model universe was not obtained. These conditional
    weights are not all-model Akaike probabilities or independent tests.
    """
    models=group["models"]
    aics=[Decimal(z["aic"]) for z in models]
    best=min(aics)
    raw=[math.exp(-float(aic-best)/2) for aic in aics]
    total=math.fsum(raw)
    out=[]
    for row,aic,rw in zip(models,aics,raw):
        delta=aic-best
        if Decimal(row["delta_aic"])!=delta:
            raise ValueError("2019 printed delta-AIC not equal to printed AIC contrast")
        out.append({
            "original_author_model":row["model"],
            "printed_AIC":row["aic"],
            "reported_delta_AIC":str(delta),
            "within_five_candidate_conditional_Akaike_weight":round(rw/total,12),
            "model_weight_is_not_gene_causality_or_heldout_prediction":True,
        })
    return {
        "response":group["response"],
        "n_listed_candidates":len(out),
        "published_best_model":"Chi",
        "published_runner_up":"F3h",
        "best_vs_runner_up_delta_AIC":str(aics[1]-aics[0]),
        "relative_best_vs_runner_up_AIC_evidence_ratio":round(raw[0]/raw[1],6),
        "relative_best_vs_runner_up_ratio_interval_for_rounded_AIC_plusminus_0p05":
            [round(math.exp((float(aics[1]-aics[0])-.1)/2),6),
             round(math.exp((float(aics[1]-aics[0])+.1)/2),6)],
        "conditional_Chi_weight_among_only_five_printed":round(raw[0]/total,12),
        "printed_AIC_values_rounded_to_one_decimal":True,
        "full_candidate_set_completeness_unverified":True,
        "model_fits_not_independently_recomputed":True,
        "all_five":out,
    }


def original_tissue_correlation(cells):
    """Independent exact two-decimal printed matrix audits, not original fit."""
    mat=[[F(x) for x in row] for row in cells]
    petal=[i for i,k in enumerate(EXPECTED_DIMENSIONS) if k.startswith("petal_")]
    leaf=[i for i,k in enumerate(EXPECTED_DIMENSIONS) if k.startswith("leaf_")]
    pairs=[(i,j,mat[i][j]) for i in range(9) for j in range(i+1,9)]
    within=[v for i,j,v in pairs if (i in petal)==(j in petal)]
    between=[v for i,j,v in pairs if (i in petal)!=(j in petal)]
    if len(within)!=18 or len(between)!=18:
        raise ArithmeticError("2016 printed 9x9 tissue-pair grouping changed")
    bykey=lambda a,b:mat[EXPECTED_DIMENSIONS.index(a)][EXPECTED_DIMENSIONS.index(b)]
    key={
        "within_petal_di_antho_vs_di_flav":bykey("petal_di_anthocyanidin","petal_di_flavonol"),
        "within_petal_tri_antho_vs_tri_flav":bykey("petal_tri_anthocyanidin","petal_tri_flavonol"),
        "within_petal_mono_flav_vs_di_flav":bykey("petal_mono_flavonol","petal_di_flavonol"),
        "within_leaf_mono_flav_vs_di_flav":bykey("leaf_mono_flavonol","leaf_di_flavonol"),
        "cross_petal_di_antho_vs_leaf_di_flav":bykey("petal_di_anthocyanidin","leaf_di_flavonol"),
        "cross_petal_di_antho_vs_leaf_mono_flav":bykey("petal_di_anthocyanidin","leaf_mono_flavonol"),
    }
    if key!={"within_petal_di_antho_vs_di_flav":F("0.68"),
            "within_petal_tri_antho_vs_tri_flav":F("0.52"),
            "within_petal_mono_flav_vs_di_flav":F("0.78"),
            "within_leaf_mono_flav_vs_di_flav":F("0.96"),
            "cross_petal_di_antho_vs_leaf_di_flav":F("0.60"),
            "cross_petal_di_antho_vs_leaf_mono_flav":F("0.54")}:
        raise ValueError("Original reported biochemical coexpression pairs drift")
    return {
        "original_published_matrix_dimensions":9,
        "published_correlation_biological_individuals":26,
        "published_correlation_species":9,
        "printed_phylogenetic_lambda":"0.342",
        "unique_within_tissue_pairs":len(within),
        "unique_between_tissue_pairs":len(between),
        "within_tissue_pairs_abs_r_at_least_0p5":sum(abs(v)>=F(1,2) for v in within),
        "between_tissue_pairs_abs_r_at_least_0p5":sum(abs(v)>=F(1,2) for v in between),
        "selected_original_printed_r":{k:str(v) for k,v in key.items()},
        "source_pairwise_significance_not_retested":True,
        "one_published_matrix_not_36_independent_plant_experiments":True,
        "flower_to_leaf_cross_tissue_r_not_zero_or_universally_strong":True,
    }


def report(source,gate,v24):
    source_verification(source,gate,v24)
    pigment={z["response"]:conditional_akaike(z)
             for z in source["source_2019"]["reported_models"]}
    corr=original_tissue_correlation(source["source_2016"][
        "original_printed_correlation_matrix"])
    if pigment["anthocyanidins"]["best_vs_runner_up_delta_AIC"]!="12.2":
        raise ArithmeticError("Original anthocyanidin Chi-F3h delta-AIC drift")
    if pigment["flavonols"]["best_vs_runner_up_delta_AIC"]!="9.9":
        raise ArithmeticError("Original flavonol Chi-F3h delta-AIC drift")
    return {
        "schema":"anto_v25_original_literature_endpoint_separation_flavonoid_evidence",
        "date":"2026-10-28",
        "real_original_published_source_tables":"2016 Berardi Table2 + 2019 Larter Table2",
        "original_2016_source":"10.1016/j.phytochem.2016.05.007",
        "original_2019_source":"10.1002/dvdy.82",
        "not_independent_plant_clades_or_same_original_individuals":True,
        "published_quantitative_outcomes":pigment,
        "published_tissue_metabolite_correlations":corr,
        "source_reported_white_fixed_vs_white_morph_gene_expression_class_contrast":{
            "original_paper_claim":"Fixed-white species show large late DFR/ANS expression reductions; most within-species acyanic morphs resemble pigmented conspecific expression profiles.",
            "source":"Larter et al. 2019 Results 3.3, Figs. 5 and 6; original author's descriptive/class-separation finding",
            "new_held_out_classification_model_fitted_by_ANTO":False,
            "does_not_equal_quantitative_pigment_PGLS_ranking":True,
        },
        "horizon_target_separation":{
            "continuous_chemical_amount":"original source-published PGLS table for anthocyanidins and flavonols; upstream CHI best among reported models",
            "categorical_fixed_white_vs_polymorphic_white_mechanism":"original authors' multivariate results; downstream DFR/ANS-related differences",
            "other_tissue_biochemistry":"2016 separately published original flower-leaf correlation; not a current full flower toolkit assay",
            "ancestral_gene_function_native_vs_donor":"NO QUALIFYING SAME-PLANT CAUSAL EXPERIMENT HERE",
            "multi_generation_trait_reexpression":"NOT OBSERVED OR ESTIMATED",
        },
        "no_original_2019_individual_HPLC_qPCR_data_replayed":True,
        "no_natural_colour_return_or_ancestral_toolkit_function_confirmed":True,
        "not_posterior_gene_causal_probability_or_cross_clade_effect":True,
        "Cirsium_organismal_empirical_owner":"zuizui0223/aza3",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path)
    p.add_argument("--expected",type=Path)
    args=p.parse_args()
    value=report(load(SOURCE),load(RAW_GATE),load(V24_TABLE))
    if args.expected and value!=load(args.expected):
        raise RuntimeError("Original published v25 AIC/correlation arithmetic changed")
    content=json.dumps(value,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding="utf-8")
    else:
        print(content)


if __name__=="__main__":
    main()
