"""Original processed-data Petunieae 2023 same-source species holdout audit.

Scientific scope:
- Wheeler et al. (2023), Proc. R. Soc. B DOI 10.1098/rspb.2023.0275,
  processed OSF zg9cu HPLC species means + normalized bud RNA expression.
- Every original file passes independently registered original OSF MD5,
  SHA256 and size before *any* CSV cell is parsed. Only derived summaries
  are emitted, not original source gene/pigment tables.
- 59 Petunieae ingroup source species; Browallia americana BROW is the
  outgroup and is excluded. 2023 author species HPLC and RNA data are
  processed means (not new raw HPLC quantification or genomic function).
- Descriptive out-of-sample species LOO and unverified-prefix-code block
  holdouts, not an independent evolutionary gain/loss lineage holdout.
- No hyperparameter tuning on the held-out response. Ridge lambda=10
  standardized log1p original TPM; intercept learned from TRAIN ONLY.
- Outcomes are log1p(sum of 6 original anthocyanidins), log1p(sum
  of 3 original flavonols), independently evaluated. Genes never include
  the target chemical columns (no direct target leakage).
- No directly measured real human-visible flower color covariate in this
  source table. Null intercept baseline is NOT a color-only model.
- Conditional single-tribe pattern, no claim of recovering ancestral
  trait after historical loss, restored inherited function, seed fitness,
  independent phylogenetic replication or macroevolutionary causal mechanism.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
from collections import defaultdict
from pathlib import Path

from audit_petunieae_2023_processed_files_v26 import (
    FILES,registered_location,verify_original,get_official_file
)

OUTGROUP="BROW"
EXPECTED_TOTAL=60
ANTHOCYANIDINS=("Del_mgg","Cyan_mgg","Pet_mgg","Pel_mgg","Peon_mgg","Malv_mgg")
FLAVONOLS=("myr_mgg","quer_mgg","kaem_mgg")
MODELS={
    "train_mean_baseline":(),
    "CHI_only":("CHI",),
    "early_three":("CHS","CHI","F3H"),
    "late_two":("DFR","ANS"),
    "branch_three":("F3primeH","F35H","FLS"),
    "joint_six":("CHI","F3primeH","F35H","DFR","ANS","FLS"),
}
RIDGE_LAMBDA=10.0


def parse_named_csv(data,expected):
    """Parse exact source CSV without treating technical rows as new plants."""
    text=data.decode("utf-8-sig")
    reader=csv.DictReader(io.StringIO(text))
    h=reader.fieldnames
    if not h or len(h)!=len(set(h)):
        raise ValueError("Original source has no rectangular, distinct column names")
    key=h[0]
    if key not in ("","species","key_0"):
        raise ValueError("Original first species identifier column unexpected")
    records={}
    for row in reader:
        if not isinstance(row,dict) or None in row:
            raise ValueError("Original source CSV not rectangular")
        code=row[key].strip()
        if not code or code in records or len(code)>12:
            raise ValueError("Original source species keys missing or duplicate")
        records[code]=row
    if len(records)!=expected:
        raise ValueError("Unexpected original 2023 source row count")
    return records,h,key


def numeric(row,key):
    val=row.get(key)
    if val is None or not val.strip():
        raise ValueError(f"Original measured predictor/response missing: {key}")
    try:
        num=float(val)
    except ValueError as exc:
        raise ValueError(f"Original processed nonnumeric {key}") from exc
    if not math.isfinite(num) or num<0:
        raise ValueError(f"Original numeric original source negative/nonfinite: {key}")
    return num


def reconcile_original(matched,expr,pigment):
    """Exact species, pigment and gene cross-source joins, preserving original keys."""
    main,hm,km=parse_named_csv(matched,EXPECTED_TOTAL)
    genes,hg,kg=parse_named_csv(expr,EXPECTED_TOTAL)
    colors,hc,kc=parse_named_csv(pigment,EXPECTED_TOTAL)
    if set(main)!=set(genes) or set(main)!=set(colors):
        raise ValueError("Original 2023 processed expression and chemical taxa do not match")
    if OUTGROUP not in main:
        raise ValueError("Original Petunieae data must keep BROW outgroup identified")
    if not set(ANTHOCYANIDINS+FLAVONOLS).issubset(hm) or not set(
        ANTHOCYANIDINS+FLAVONOLS).issubset(hc):
        raise ValueError("Original HPLC species chemical endpoints missing")
    source_genes=set(hg)&set(hm)-set(ANTHOCYANIDINS+FLAVONOLS)
    if not set().union(*map(set,MODELS.values())).issubset(source_genes):
        raise ValueError("Predeclared molecular predictors not in both source tables")
    source_pigment_cells=0
    source_gene_cells=0
    for code in main:
        for pigment_key in ANTHOCYANIDINS+FLAVONOLS:
            a=numeric(main[code],pigment_key)
            b=numeric(colors[code],pigment_key)
            if not math.isclose(a,b,abs_tol=1e-10,rel_tol=1e-10):
                raise ValueError("Original same-species HPLC rows disagree between author source files")
            source_pigment_cells+=1
        for name in set().union(*map(set,MODELS.values())):
            a=numeric(main[code],name)
            b=numeric(genes[code],name)
            if not math.isclose(a,b,abs_tol=1e-7,rel_tol=1e-9):
                raise ValueError("Original same-species gene rows disagree between author source files")
            source_gene_cells+=1
    ingroup=sorted(set(main)-{OUTGROUP})
    if len(ingroup)!=59:
        raise ValueError("Original source 59 Petunieae plus outgroup denominator changed")
    samples={}
    for code in ingroup:
        x={name:math.log1p(numeric(main[code],name))
           for name in set().union(*map(set,MODELS.values()))}
        responses={
            "anthocyanidins_log1p_total":math.log1p(
                math.fsum(numeric(main[code],x) for x in ANTHOCYANIDINS)),
            "flavonols_log1p_total":math.log1p(
                math.fsum(numeric(main[code],x) for x in FLAVONOLS)),
        }
        samples[code]={"genes":x,"responses":responses}
    return samples,{
        "original_author_processed_species_rows":len(main),
        "Petunieae_ingroup_species_rows":len(ingroup),
        "Browallia_americana_outgroup_removed":True,
        "source_species_set_match_all_three_files":True,
        "source_original_pigment_cells_verified_by_independent_table":source_pigment_cells,
        "source_original_gene_cells_verified_by_independent_table":source_gene_cells,
        "molecular_predictors_prespecified":sorted(set().union(*map(set,MODELS.values()))),
        "no_visible_hue_observed_baseline_in_this_table":True,
    }


def linear_solve(a,b):
    """Independent small dense ridge system via deterministic partial pivot."""
    n=len(b)
    if n==0:
        return []
    m=[list(a[i])+[b[i]] for i in range(n)]
    for j in range(n):
        pivot=max(range(j,n),key=lambda k:abs(m[k][j]))
        if abs(m[pivot][j])<1e-12:
            raise ValueError("Singular ridge system, result cannot be admitted")
        m[j],m[pivot]=m[pivot],m[j]
        v=m[j][j]
        for c in range(j,n+1):
            m[j][c]/=v
        for i in range(n):
            if i==j:continue
            fac=m[i][j]
            for c in range(j,n+1):
                m[i][c]-=fac*m[j][c]
    return [m[i][n] for i in range(n)]


def ridge_predict(train,test,endpoint,predictors,lam=RIDGE_LAMBDA):
    """All training means and scales are fit *within each holdout fold*."""
    if len(train)<3:
        raise ValueError("At least three training species required")
    if not predictors:
        return math.fsum(train[k]["responses"][endpoint] for k in train)/len(train)
    n=len(train)
    names=tuple(predictors)
    y=[train[k]["responses"][endpoint] for k in train]
    ymean=math.fsum(y)/n
    xlog=[[train[k]["genes"][name] for name in names] for k in train]
    means=[math.fsum(row[j] for row in xlog)/n for j in range(len(names))]
    sd=[math.sqrt(math.fsum((row[j]-means[j])**2 for row in xlog)/n)
        for j in range(len(names))]
    scale=[s if s>1e-10 else 1.0 for s in sd]
    x=[[(row[j]-means[j])/scale[j] for j in range(len(names))] for row in xlog]
    kdim=len(names)
    gram=[[math.fsum(row[i]*row[j] for row in x)+(lam if i==j else 0)
           for j in range(kdim)] for i in range(kdim)]
    right=[math.fsum(row[j]*(yi-ymean) for row,yi in zip(x,y))
           for j in range(kdim)]
    beta=linear_solve(gram,right)
    zs=[(test["genes"][name]-means[j])/scale[j] for j,name in enumerate(names)]
    return ymean+math.fsum(b*z for b,z in zip(beta,zs))


def blocked_folds(samples,kind):
    ids=sorted(samples)
    if kind=="species_LOO":
        return [[code] for code in ids]
    if kind=="first_two_source_code_characters_NOT_accepted_genus":
        groups=defaultdict(list)
        for code in ids:
            groups[code[:2]].append(code)
        if len(groups)<3 or max(map(len,groups.values()))>len(ids)-5:
            raise ValueError("Source identifier prefix grouping has insufficient independent training")
        return [groups[k] for k in sorted(groups)]
    raise ValueError("Unregistered validation fold type")


def evaluate(samples,kind):
    folds=blocked_folds(samples,kind)
    endpoints=("anthocyanidins_log1p_total","flavonols_log1p_total")
    predictions={t:{m:{} for m in MODELS} for t in endpoints}
    ids=set(samples)
    for holdout in folds:
        held=set(holdout)
        train=sorted(ids-held)
        if len(train)<5:
            raise ValueError("Insufficient training species for original source holdout")
        for endpoint in endpoints:
            for name,genes in MODELS.items():
                for code in holdout:
                    predictions[endpoint][name][code]=ridge_predict(
                        {k:samples[k] for k in train},samples[code],endpoint,genes)
    out={}
    for endpoint in endpoints:
        truth={x:samples[x]["responses"][endpoint] for x in sorted(ids)}
        baseline=None
        metrics={}
        n=len(ids)
        variance=math.fsum((y-math.fsum(truth.values())/n)**2 for y in truth.values())
        for name in MODELS:
            if set(predictions[endpoint][name])!=ids:
                raise ArithmeticError("Each real source species must have one out-of-fold prediction")
            errors=[truth[code]-predictions[endpoint][name][code]
                    for code in sorted(ids)]
            se=math.fsum(e*e for e in errors)
            ae=math.fsum(abs(e) for e in errors)
            rmse=math.sqrt(se/n)
            if name=="train_mean_baseline":
                baseline=rmse
            metrics[name]={
                "n_tested_original_ingroup_species":n,
                "out_of_fold_RMSE_log1p_mg_per_g":round(rmse,8),
                "out_of_fold_MAE_log1p_mg_per_g":round(ae/n,8),
                "out_of_fold_R2_vs_global_mean_not_phylogenetic":round(1-se/variance,8)
                  if variance>0 else None,
                "relative_RMSE_to_source_LOO_baseline":None if baseline is None
                   else round(rmse/baseline,8) if baseline>0 else None,
                "strict_test_species_outside_own_train_set":True,
            }
        out[endpoint]={
           "source_original_zero_total_species":sum(
             math.expm1(z)==0 for z in truth.values()),
           "phenotype_mean_log1p_mg_per_g":round(math.fsum(truth.values())/n,8),
           "pigment_chemical_response_is_NOT_visual_hue_or_ancestral_reexpression":True,
           "fixed_unretuned_ridge_lambda":RIDGE_LAMBDA,
           "models":metrics,
        }
    return {
        "holdout_type":kind,
        "n_original_ingroup_species":len(ids),
        "n_folds":len(folds),
        "fold_test_sizes_sorted":sorted(map(len,folds)),
        "code_prefix_is_NOT_verified_taxonomic_genus_or_an_independent_regain_clade":True,
        "refitted_train_only_scaling_and_intercept":True,
        "endpoints":out,
    }


def source_real_data(fetch=get_official_file):
    needed=("expression_pigment_matched","expression_gene_mean","pigment_mean")
    original={}
    for key in needed:
        row=FILES[key]
        data=fetch(registered_location(row))
        verify_original(data,row)
        original[key]=data
    return reconcile_original(original["expression_pigment_matched"],
                              original["expression_gene_mean"],
                              original["pigment_mean"])


def report(fetch=get_official_file):
    samples,origin=source_real_data(fetch)
    validations={
        kind:evaluate(samples,kind)
        for kind in ("species_LOO","first_two_source_code_characters_NOT_accepted_genus")
    }
    return {
        "schema":"anto_v26_original_2023_petunieae_species_pigment_expression_prediction",
        "source_study_DOI":"10.1098/rspb.2023.0275",
        "original_OSF_project":"zg9cu",
        "original_author_processed_flower_RNA_and_HPLC":True,
        "source_file_hashes_verified":{
            name:{"source_registered_md5":FILES[name]["md5"],
                  "source_registered_sha256":FILES[name]["sha256"],
                  "source_declared_bytes":FILES[name]["size"]}
            for name in ("expression_pigment_matched","expression_gene_mean","pigment_mean")
        },
        "source_reconciled":origin,
        "within_one_original_clade_and_author_dataset":True,
        "original_sra_RNA_biological_requantification_not_performed":True,
        "not_a_phenotypic_visible_color_baseline_or_future_regain":True,
        "fixed_models_no_response_specific_feature_search_or_tuning":list(MODELS),
        "no_phylogenetically_independent_real_ancestral_recovery":True,
        "no_full_genetic_competence_or_causal_intervention_measured":True,
        "no_external_clade_transfer_or_transgenerational_fitness_validation":True,
        "Cirsium_organismal_empirical_owner":"zuizui0223/aza3",
        "validation":validations,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--expected",type=Path)
    args=p.parse_args()
    result=report()
    if args.expected is not None:
        reference=json.loads(args.expected.read_text(encoding="utf-8"))
        if result!=reference:
            raise RuntimeError("Original 2023 processed Petunieae out-of-fold performance changed")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ANTO_V26_ORIGINAL_SPECIES_HOLDOUT_SUMMARY",json.dumps({
        "original_species":result["source_reconciled"]["Petunieae_ingroup_species_rows"],
        "original_source_genes_reconciled":result["source_reconciled"]["source_original_gene_cells_verified_by_independent_table"],
        "original_source_pigment_cells_reconciled":result["source_reconciled"]["source_original_pigment_cells_verified_by_independent_table"],
        "endpoints":{typ:{
             response:{"models":{k:{"rmse":v["out_of_fold_RMSE_log1p_mg_per_g"],
                                    "rmse_over_null":v["relative_RMSE_to_source_LOO_baseline"]}
                               for k,v in info["models"].items()}}
             for response,info in check["endpoints"].items()}
           for typ,check in result["validation"].items()}
    },sort_keys=True))


if __name__=="__main__":
    main()
