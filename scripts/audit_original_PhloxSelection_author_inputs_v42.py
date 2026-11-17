"""ANTO v42: checksum-pinned original AUTHOR code input replay, not full MCMC.

Read PhloxHopkins/PhloxSelection 2014 commit original gene-flow simulation
*example input* at a user-provided checked-out path. Verifies original Git
blob hashes for source, 32 population phenotype vectors, family-by-maternal
phenotype source input, and provenance. Prints aggregate-only derived values.
Original raw author source files never redistributed by ANTO.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/"data/external/anto_v42_original_PhloxSelection_author_git_manifest_20261117.json"
V41=ROOT/"data/external/anto_v41_Phlox2014_original_spatial_study_tableS1_provisional_20261117.json"
COLORS=("DB","DR","LB","LR")
MOTHER_LABEL={0:"DB",1:"DR",2:"LB",3:"LR"}
AUTHOR_COMMIT="60ff508170f1be592bbc6fe7ed46ac6e5b459026"
AUTHOR_REPO="PhloxHopkins/PhloxSelection"

class OriginalSourceError(ValueError):
    pass


def git_blob_sha1(payload):
    if not isinstance(payload,bytes):
        raise OriginalSourceError("GIT_BLOB_REQUIRES_ORIGINAL_BYTES")
    return hashlib.sha1(b"blob "+str(len(payload)).encode("ascii")+b"\0"+payload).hexdigest()


def require(cond,message):
    if not cond:
        raise OriginalSourceError(message)


def digest_file(path,sha):
    payload=path.read_bytes()
    actual=git_blob_sha1(payload)
    require(actual==sha,"AUTHOR_ORIGINAL_GIT_BLOB_SHA1_DRIFT_"+path.name)
    return payload


def _parse_table(raw,expected_header):
    try:
        text=raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise OriginalSourceError("AUTHOR_INPUT_NOT_UTF8") from exc
    reader=csv.reader(io.StringIO(text),delimiter="\t")
    header=next(reader,None)
    require(header==list(expected_header),"AUTHOR_INPUT_HEADER_SCHEMA_DRIFT")
    result=[]
    for line,row in enumerate(reader,2):
        if not row or all(not c.strip() for c in row):
            continue # author Fam_T4 includes trailing fully blank padding lines
        require(len(row)==len(expected_header) and all(c.strip() for c in row),
                "AUTHOR_INPUT_NONEMPTY_BAD_COLUMN_COUNT_LINE_"+str(line))
        result.append([x.strip() for x in row])
    require(len(result)>0,"AUTHOR_INPUT_ZERO_VALID_ROWS")
    return result


def nat(x):
    if not x.isdigit() or (len(x)>1 and x.startswith("+")):
        raise OriginalSourceError("AUTHOR_INPUT_NONINTEGER_OR_NEGATIVE_COUNT")
    return int(x)


def parse_phen(raw):
    rows=_parse_table(raw,("pop","distance","DB","DR","LB","LR"))
    out=[]
    seen=set()
    previous_d=-math.inf
    for row in rows:
        pop=nat(row[0])
        require(pop not in seen,"AUTHOR_PHEN_DUPLICATE_POPULATION_ID")
        seen.add(pop)
        try:
            distance=float(row[1])
        except ValueError as exc:
            raise OriginalSourceError("AUTHOR_PHEN_NONNUMERIC_DISTANCE") from exc
        require(math.isfinite(distance) and distance>=0 and distance>=previous_d,
                "AUTHOR_PHEN_INVALID_OR_UNSORTED_DISTANCE")
        previous_d=distance
        nums=[nat(v) for v in row[2:]]
        require(sum(nums)>0,"AUTHOR_PHEN_SOURCE_ZERO_SAMPLE_POP")
        out.append({"population":pop,"distance_km":distance,
                    "by_color":dict(zip(COLORS,nums))})
    return out


def parse_fam(raw,valid_pops):
    rows=_parse_table(raw,("population","color","DB","DR","LB","LR"))
    out=[]
    for row in rows:
        pop=nat(row[0]);mother=nat(row[1])
        require(pop in valid_pops,"AUTHOR_FAM_POP_NOT_IN_SAME_TRANSECT")
        require(mother in MOTHER_LABEL,"AUTHOR_MOTHER_COLOR_CODE_CHANGED")
        nums=[nat(x) for x in row[2:]]
        require(sum(nums)>0,"AUTHOR_FAM_ZERO_CHILDREN_NOT_CONFIRMED_MOTHER")
        out.append({"population":pop,
                    "maternal_phenotype":MOTHER_LABEL[mother],
                    "offspring_by_color":dict(zip(COLORS,nums)),
                    "n_offspring":sum(nums)})
    return out


def parent_child_matrix(families):
    return {
        c:{
            "n_source_family_rows":sum(f["maternal_phenotype"]==c for f in families),
            "offspring_n":sum(f["n_offspring"] for f in families if f["maternal_phenotype"]==c),
            "children_by_colour":{x:sum(f["offspring_by_color"][x] for f in families
                                     if f["maternal_phenotype"]==c)
                                  for x in COLORS},
        }
        for c in COLORS
    }


def report(original_checkout,manifest,v41):
    require(manifest.get("schema")==
            "anto_v42_source_locked_original_author_PhloxSelection_input_inventory",
            "V42_SOURCE_MANIFEST_INVALID")
    require(manifest.get("canonical_GitHub_repo")==AUTHOR_REPO and
            manifest.get("original_commit_sha")==AUTHOR_COMMIT and
            manifest.get("original_author_README_calls_files_example_input") is True,
            "V42_ORIGINAL_AUTHOR_SOURCE_IDENTITY_DRIFT")
    require(v41.get("schema")=="anto_v41_original2014_Phlox_spatial_population_selection_source_registry",
            "V41_PREVIOUS_DIGITIZED_SOURCE_SCHEMA_DRIFT")
    try:
        head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=original_checkout,
                        stderr=subprocess.PIPE,text=True).strip()
    except (OSError,subprocess.CalledProcessError) as exc:
        raise OriginalSourceError("NOT_A_GIT_CHECKOUT_OF_AUTHOR_SOURCE") from exc
    require(head==AUTHOR_COMMIT,"AUTHOR_SOURCE_GIT_COMMIT_NOT_EXACT_PIN")

    blobs=manifest["all_author_files_fixed_git_blob_sha1"]
    expected_files=(
        ["README.md","code/WorldOutput.cpp","code/main.cpp",
         "parameterFiles/inPhlox_nuissance.params","parameterFiles/seedSet.params"]
        +[f"inputFiles/inPhlox_Phen_T{i}.data" for i in range(1,6)]
        +[f"inputFiles/inPhlox_Fam_T{i}.data" for i in range(1,6)]
        +[f"inputFiles/init_Phlox_IF{i}.txt" for i in range(1,6)]
    )
    require(set(blobs)==set(expected_files),"AUTHOR_SOURCE_REGISTERED_FILE_SET_DRIFT")
    data={path:digest_file(original_checkout/path,blobs[path]) for path in expected_files}
    readme=data["README.md"].decode("utf-8-sig")
    require("included example input files" in readme and
            'observed offspring counts from each maternal plant' in readme and
            "inPhlox_Phen_T" in readme and "inPhlox_Fam_T" in readme,
            "ORIGINAL_README_EXAMPLE_SOURCE_GRANULARITY_LOST")
    main=data["code/main.cpp"].decode("utf-8-sig")
    require("#define mcmcRealz 75000" in main,
            "ORIGINAL_MAIN_LONG_MCMC_RUNTIME_GATE_DRIFT")
    code=data["code/WorldOutput.cpp"].decode("utf-8-sig")
    require("Phenotype freq: DB " in code and
            "pheno[0]" in code and "pheno[1]" in code and
            "pheno[2]" in code and "pheno[3]" in code,
            "MOTHER_SOURCE_PHENOTYPE_NUMERIC_INDEX_UNKNOWN")
    require(manifest["original_code_maternal_coding_0DB_1DR_2LB_3LR"] is True,
            "MOTHER_COLOR_CODE_REGISTER_DRIFT")

    expected_v41=v41["source_transcription"]["population_phenotype_counts"]
    require(len(expected_v41)==32,"V41_EXPECTED_32_POPULATIONS_MISSING")
    fen=[]
    fam=[]
    per=[]
    for t in range(1,6):
        ph=f"inputFiles/inPhlox_Phen_T{t}.data"
        ff=f"inputFiles/inPhlox_Fam_T{t}.data"
        init=f"inputFiles/init_Phlox_IF{t}.txt"
        require(data[ph]==data[init],
                "AUTHOR_INIT_AND_PHEN_NOT_DUPLICATED_SOURCE_"+str(t))
        p=parse_phen(data[ph])
        valid={x["population"] for x in p}
        f=parse_fam(data[ff],valid)
        fen.extend((t,x) for x in p)
        fam.extend((t,x) for x in f)
        adult=sum(sum(x["by_color"].values()) for x in p)
        kids=sum(x["n_offspring"] for x in f)
        per.append({
            "transect":t,
            "n_adult_populations":len(p),
            "n_original_family_rows":len(f),
            "n_offspring_counted":kids,
            "n_adult_phenotype_counts":adult,
            "family_source_populations":len({r["population"] for r in f}),
            "phen_source_first_and_last_populations":[p[0]["population"],p[-1]["population"]],
        })
    require(len(fen)==32 and len({x["population"] for _,x in fen})==32,
            "AUTHOR_CODE_MUST_HAVE_32_UNIQUE_POPULATIONS")
    reconstructed=[]
    for t,p in fen:
        reconstructed.append([t,p["population"],p["by_color"]["LB"],
                              p["by_color"]["DB"],p["by_color"]["LR"],p["by_color"]["DR"]])
    require(reconstructed==expected_v41,
            "ORIGINAL_AUTHOR_S1_AND_V41_MANUAL_TABLE_ANY_ROW_MISMATCH")
    adult_by_color={c:sum(x["by_color"][c] for _,x in fen) for c in COLORS}
    mother_mat=parent_child_matrix([x for _,x in fam])
    offspring_by_color={c:sum(x["offspring_by_color"][c] for _,x in fam) for c in COLORS}
    adult_total=sum(adult_by_color.values())
    child_total=sum(offspring_by_color.values())
    family_total=len(fam)
    require(adult_total==11379 and len(fam)==371 and child_total==2496,
            "2014_AUTHOR_EXAMPLE_FULL_COHORT_COUNTER_DRIFT")
    require(adult_by_color=={"DB":1396,"DR":4607,"LB":4292,"LR":1084},
            "AUTHOR_S1_CLASS_COUNTS_SOURCE_DRIFT")
    require(offspring_by_color=={"DB":375,"DR":776,"LB":1093,"LR":252},
            "AUTHOR_FAM_OFFSPRING_PUBLISHED_CODE_COUNTS_DRIFT")
    require({c:mother_mat[c]["n_source_family_rows"] for c in COLORS}==
            {"DB":59,"DR":128,"LB":151,"LR":33},
            "AUTHOR_FAM_MATERNAL_COLOR_COUNTS_DRIFT")
    adult_pop_ids={x["population"] for _,x in fen}
    mother_pop_ids={x["population"] for _,x in fam}
    missing=sorted(adult_pop_ids-mother_pop_ids)
    require(missing==[330,332],"AUTHOR_EXAMPLE_MATERNAL_POP_MISSING_SOURCE_DRIFT")
    expected=manifest["original_article_claimed_counts"]
    require(expected=={"adult":12024,"maternal_plants":382,
                       "offspring":2581,"populations":32,"transects":5},
            "PUBLISHED_MAIN_PAPER_OBSERVATION_DENOMINATOR_DRIFT")
    missing_counts={
        "adult_main_minus_example":expected["adult"]-adult_total,
        "maternal_main_minus_example":expected["maternal_plants"]-family_total,
        "offspring_main_minus_example":expected["offspring"]-child_total,
    }
    require(missing_counts=={"adult_main_minus_example":645,
                             "maternal_main_minus_example":11,
                             "offspring_main_minus_example":85},
            "ORIGINAL_PAPER_AUTHOR_EXAMPLE_UNRECONCILED_DIFFERENCES_DRIFT")

    seedline=data["parameterFiles/seedSet.params"].decode("utf-8").strip().split()
    require(len(seedline)==13,"ORIGINAL_START_PARAMETERS_LENGTH_DRIFT")
    start=[float(x) for x in seedline]
    require(all(math.isfinite(v) for v in start),"ORIGINAL_START_PARAMETERS_NONFINITE")
    require(start==[.279,.362,.612,.639,.362,0,.72,.88,.91,.90,3.7,6.3,5.5],
            "SOURCE_2014_AUTHOR_INITIAL_MCMC_NOT_ORIGINAL_FITTED_TABLE1")

    return {
        "schema":"anto_v42_original_author_PhloxSelection_example_inputs_integrity_and_population_offspring_counts",
        "original_author_source":{
            "publication_DOI":manifest["original_study_DOI"],
            "GitHub_repository":AUTHOR_REPO,
            "source_supplement_links_GitHub_at_historical_owner":True,
            "original_code_commit_pinned":head,
            "original_git_blob_sha1_verification_count":len(expected_files),
            "all_source_digests_verified_from_pinned_original_checkout":True,
            "original_author_readme_labels_inputs_as_examples":True,
            "main_source_default_75000_MCMC_steps_not_executed":True,
            "phen_and_init_same_five_files_not_independent_data":True,
            "original_source_code_not_copied_to_ANTO":True,
        },
        "source_2014_all_32_S1_phenotype_rows_exactly_match_independent_v41_manual_transcription":True,
        "n_adult_populations":len(fen),
        "adult_phenotype_total_original_author_code_example":adult_total,
        "adult_by_color":adult_by_color,
        "mother_source_family_row_total":family_total,
        "mother_source_rows_without_unique_maternal_plant_ID":True,
        "source_offspring_phenotype_total":child_total,
        "source_offspring_by_color":offspring_by_color,
        "source_mother_by_offspring_color":mother_mat,
        "per_transect_source_counts":per,
        "n_source_populations_with_mothers":len(mother_pop_ids),
        "source_populations_with_no_example_maternal_records":missing,
        "source_adult_populations_that_have_no_mother_row_are_not_absent_from_manuscript":True,
        "independent_ancestor_gene_phenotype_fitness_causality_newly_measured":False,
        "main_article_reported_original_observation_counts":expected,
        "published_main_minus_author_example_input":missing_counts,
        "source_cohort_reconciliation_state":
            "AUTHOR_INPUT_EXAMPLE_DOES_NOT_RECONCILE_PUBLISHED_MAIN_COHORT",
        "entire_original_manuscript_2014_MCMC_likelihood_reproduced":False,
        "maternal_phenotype_not_original_maternal_genotype_haplotype":True,
        "sibling_offspring_not_independent_evolutionary_lineages":True,
        "source_mother_codes":{str(k):v for k,v in MOTHER_LABEL.items()},
        "original_2014_initial_parameter_values_not_fitted_original_posterior":start,
        "no_arbitrary_filling_of_645_11_85_missing_observation_counts":True,
        "no_new_Phylo_Phlox_allele_selection_or_regain_probability":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--author-git-checkout",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--expected",type=Path)
    a=p.parse_args()
    obj=report(a.author_git_checkout,
               json.loads(MANIFEST.read_text(encoding="utf-8")),
               json.loads(V41.read_text(encoding="utf-8")))
    if a.expected and json.loads(a.expected.read_text(encoding="utf-8"))!=obj:
        raise OriginalSourceError("V42_SOURCE_PERSISTENT_DERIVED_RESULT_DRIFT")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(obj,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("ANTO_V42_ORIGINAL_AUTHOR_INPUT_SOURCE_REPLAY",json.dumps({
        "original_commit":obj["original_author_source"]["original_code_commit_pinned"],
        "all_S1_population_colour_rows_v41_exact":obj[
            "source_2014_all_32_S1_phenotype_rows_exactly_match_independent_v41_manual_transcription"],
        "author_input_adults":obj["adult_phenotype_total_original_author_code_example"],
        "author_input_mothers":obj["mother_source_family_row_total"],
        "author_input_children":obj["source_offspring_phenotype_total"],
        "main_minus_example":obj["published_main_minus_author_example_input"],
        "missing_mother_pop":obj["source_populations_with_no_example_maternal_records"],
        "biological_full_original_mcmc_refit":False,
    },sort_keys=True))
    packed=json.dumps(obj,sort_keys=True,separators=(",",":"))
    block=6000;n=(len(packed)+block-1)//block
    for i in range(n):
        print("ANTO_V42_DERIVED_CHUNK",i+1,n,packed[i*block:(i+1)*block])


if __name__=="__main__":
    main()
