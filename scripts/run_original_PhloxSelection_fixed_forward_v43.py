"""Compile original unchanged PhloxSelection 2014 forward C++ and evaluate once.

The original Git tree is SHA-pinned, author's data example cohort is validated
with v42 exact original source. Only ANTO-authored temp type aliases + driver,
plus removal of all-blank padded family lines in a temporary working copy.
No MCMC, original C++ code rewriting or raw source redistribution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import audit_original_PhloxSelection_author_inputs_v42 as prev

ROOT=Path(__file__).resolve().parents[1]
DRIVER=ROOT/"scripts/PhloxFixedSourceForwardDriver_v43.cpp"
V42_RESULT=ROOT/"results/anto_v42_original_author_input_replay_20261117.json"
ORIGINAL_CPP=(
    "World.cpp","Deme.cpp","initialfreq.cpp","simulate.cpp",
    "WorldOutput.cpp","Params.cpp","utilities.cpp","likelihood.cpp",
)
COMPAT_HEADER=(
    "#ifndef ANTO_V43_MINIMAL_NR3_TYPE_ALIAS\n"
    "#define ANTO_V43_MINIMAL_NR3_TYPE_ALIAS\n"
    "#include <vector>\n"
    "typedef std::vector<double> VecDoub;\n"
    "typedef const VecDoub& VecDoub_I;\n"
    "typedef double Doub;\n"
    "#endif\n"
)
EXPECTED_SOURCE_START=[.279,.362,.612,.639,.362,0.,.72,.88,.91,.90,3.7,6.3,5.5]
EXPECTED_WORLD_COUNTS={
    1:(9,4562,121,649),2:(7,1744,87,630),3:(6,2300,56,407),
    4:(6,1424,78,586),5:(4,1349,29,224),
}
TAG_WORLD="ANTO_V43_ORIGINAL_SOURCE_TRANSECT"
TAG_TOTAL="ANTO_V43_ORIGINAL_SOURCE_TOTAL"


class FixedForwardError(ValueError):
    pass


def require(ok,what):
    if not ok:
        raise FixedForwardError(what)


def strip_only_blank_original_family_padding(raw):
    """Do not alter any nonblank original source row bytes."""
    lines=raw.splitlines(keepends=True)
    kept=[line for line in lines if line.strip()]
    rejected=[line for line in lines if not line.strip()]
    require(kept and kept[0].startswith(b"population\tcolor\tDB\tDR\tLB\tLR"),
            "AUTHOR_FAMILY_HEADER_NOT_ADMITTED")
    require(all(not x.strip() for x in rejected),
            "NOT_ONLY_BLANK_PADDING_REMOVED")
    merged=b"".join(kept)
    require(len(merged)<=len(raw),"SANITIZED_COPY_SIZE_GAIN")
    return merged,len(rejected)


def parse_original_driver_output(stdout,v42):
    """Reject partial source result or invalid probability/likelihood records."""
    world={}
    total=None
    for line in stdout.splitlines():
        if line.startswith(TAG_WORLD+" "):
            parts=line[len(TAG_WORLD)+1:].split()
            require(len(parts)==18,"DRIVER_OUTPUT_TRANSECT_FIELD_COUNT_DRIFT")
            try:
                t,pops,adults,moms,children,gens=(int(v) for v in parts[:6])
                ll,la,lf,error=(float(v) for v in parts[6:10])
                pred_adult=[float(x) for x in parts[10:14]]
                pred_child=[float(x) for x in parts[14:18]]
            except (ValueError,OverflowError) as exc:
                raise FixedForwardError("DRIVER_OUTPUT_NONNUMERIC") from exc
            require(t in EXPECTED_WORLD_COUNTS and t not in world,
                    "DRIVER_MISSING_OR_DUPLICATE_TRANSECT")
            require((pops,adults,moms,children)==EXPECTED_WORLD_COUNTS[t],
                    "SOURCE_TRANSECT_POP_MOTHER_OFFSPRING_DENOMINATOR_DRIFT")
            require(gens>0 and gens<100000,"NO_ORIGINAL_MODEL_EQUILIBRIUM_CONFIRMATION")
            require(all(math.isfinite(x) for x in [ll,la,lf,error]+pred_adult+pred_child)
                    and ll<=0 and la<=0 and lf<=0 and error<=1e-6 and error>=0,
                    "DRIVER_INVALID_OR_UNNORMALIZED_PREDICTION")
            require(math.isclose(ll,la+lf,rel_tol=0,abs_tol=1e-6),
                    "ORIGINAL_PROGRAM_LIKELIHOOD_NOT_MATCHED_BY_COMPONENTS")
            for probs in (pred_adult,pred_child):
                require(all(-1e-10<=x<=1+1e-10 for x in probs)
                        and math.isclose(sum(probs),1,rel_tol=0,abs_tol=1e-6),
                        "AGGREGATED_ORIGINAL_PREDICTED_FREQUENCIES_NOT_NORMALIZED")
            world[t]={
                "transect":t,"n_populations":pops,"n_adult_phenotype_observations":adults,
                "n_maternal_family_rows":moms,"n_offspring_phenotypes":children,
                "source_model_generations_until_convergence":gens,
                "original_model_log_likelihood":ll,
                "independent_origin_program_flower_component_log_like":la,
                "independent_origin_program_family_component_log_like":lf,
                "maximum_component_probability_row_sum_deviation":error,
                "observed_adult_number_weighted_expected_phenotype_fraction_DB_DR_LB_LR":pred_adult,
                "observed_child_number_weighted_expected_progeny_fraction_DB_DR_LB_LR":pred_child,
            }
        elif line.startswith(TAG_TOTAL+" "):
            parts=line[len(TAG_TOTAL)+1:].split()
            require(len(parts)==8 and total is None,
                    "DRIVER_TOTAL_RECORD_MISSING_OR_DUPLICATED")
            try:
                cnt=[int(x) for x in parts[:4]]
                nums=[float(x) for x in parts[4:]]
            except (ValueError,OverflowError) as exc:
                raise FixedForwardError("DRIVER_TOTAL_RECORD_NOT_NUMERIC") from exc
            require(cnt==[32,11379,371,2496],
                    "ORIGINAL_SOURCE_FULL_EXAMPLE_INPUT_DENOMINATOR_DRIFT")
            total={
                "n_populations":cnt[0],"n_adults":cnt[1],
                "n_mothers":cnt[2],"n_offspring":cnt[3],
                "original_total_log_likelihood":nums[0],
                "original_model_adult_component_log_likelihood":nums[1],
                "original_model_family_component_log_likelihood":nums[2],
                "max_normalization_error":nums[3],
            }
    require(set(world)==set(range(1,6)) and total is not None,
            "SOURCE_ORIGINAL_MODEL_NOT_ALL_FIVE_TRANSECTS_EVALUATED")
    for name,key in [
        ("n_populations","n_populations"),("n_adults","n_adult_phenotype_observations"),
        ("n_mothers","n_maternal_family_rows"),("n_offspring","n_offspring_phenotypes")]:
        require(total[name]==sum(x[key] for x in world.values()),
                "ORIGINAL_SOURCE_AGGREGATE_DENOMINATOR_DRIFT")
    for total_key,part_key in [
        ("original_total_log_likelihood","original_model_log_likelihood"),
        ("original_model_adult_component_log_likelihood",
         "independent_origin_program_flower_component_log_like"),
        ("original_model_family_component_log_likelihood",
         "independent_origin_program_family_component_log_like")]:
        require(math.isclose(total[total_key],
                             sum(x[part_key] for x in world.values()),
                             abs_tol=1e-6,rel_tol=0),
                "ORIGINAL_TOTAL_SOURCE_LIKELIHOOD_COMPONENTS_DISAGREE")
    require(math.isclose(total["original_total_log_likelihood"],
                         total["original_model_adult_component_log_likelihood"]+
                         total["original_model_family_component_log_likelihood"],
                         rel_tol=0,abs_tol=1e-6),
            "ORIGINAL_FIT_COMPONENT_SUM_DISAGREEMENT")
    require(total["max_normalization_error"]<=1e-6 and
            total["max_normalization_error"]>=0,
            "ORIGINAL_SOURCE_PROBABILITY_ROW_GATE_FAILED")
    require(v42["adult_phenotype_total_original_author_code_example"]==11379
            and v42["mother_source_family_row_total"]==371
            and v42["source_offspring_phenotype_total"]==2496,
            "V42_PAST_FROZEN_AUTHOR_COHORT_CHANGED")
    return {"per_transect":[world[t] for t in sorted(world)],"all_transects":total}


def _run_cmd(args,cwd,timeout=140):
    try:
        p=subprocess.run(args,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                         text=True,timeout=timeout,check=False)
    except subprocess.TimeoutExpired as exc:
        raise FixedForwardError("ORIGINAL_CXX_FORWARD_TECHNICAL_TIMEOUT") from exc
    if p.returncode!=0:
        raise FixedForwardError(
            "ORIGINAL_CXX_FORWARD_TECHNICAL_FAILURE_EXIT_"+str(p.returncode)+
            " STDERR="+p.stderr[-1400:].strip()+
            " STDOUT_END="+p.stdout[-900:].strip()
        )
    return p


def report(author_git_path,compiler="g++",timeout_seconds=140):
    manifest=json.loads(prev.MANIFEST.read_text())
    v41=json.loads(prev.V41.read_text())
    v42_existing=json.loads(V42_RESULT.read_text())
    replay=prev.report(author_git_path,manifest,v41)
    require(replay==v42_existing,
            "V42_ORIGINAL_SOURCE_OR_FROZEN_PRIOR_RESULT_DRIFT")
    require(not (author_git_path/"code/nr3.h").exists()
            and not (author_git_path/"nr3.h").exists()
            and not (author_git_path/"Rseeds.txt").exists(),
            "OLD_SOURCE_DEPENDENCY_MISSING_STATUS_CHANGED")
    require("Mchain" not in DRIVER.read_text() and
            "75000" not in DRIVER.read_text(),
            "ORIGINAL_MCMC_75000_STEPS_MUST_NOT_BE_IN_DRIVER")
    with tempfile.TemporaryDirectory(prefix="anto_v43_fixed_cpp_") as base:
        base=Path(base)
        work=base/"work"
        shim=base/"shim"
        work.mkdir();shim.mkdir()
        (shim/"nr3.h").write_text(COMPAT_HEADER)
        blank_removed={}
        for t in range(1,6):
            for suffix in ("Phen","Fam"):
                name=f"inPhlox_{suffix}_T{t}.data"
                path=author_git_path/"inputFiles"/name
                raw=prev.digest_file(path,manifest[
                    "all_author_files_fixed_git_blob_sha1"][f"inputFiles/{name}"])
                if suffix=="Fam":
                    content,nblank=strip_only_blank_original_family_padding(raw)
                    blank_removed[name]=nblank
                else:
                    content=raw
                (work/name).write_bytes(content)
            init=f"init_Phlox_IF{t}.txt"
            raw=prev.digest_file(author_git_path/"inputFiles"/init,
                                 manifest["all_author_files_fixed_git_blob_sha1"][
                                     f"inputFiles/{init}"])
            (work/init).write_bytes(raw)
        for name in ("inPhlox_nuissance.params","seedSet.params"):
            rel="parameterFiles/"+name
            raw=prev.digest_file(author_git_path/rel,
                                 manifest["all_author_files_fixed_git_blob_sha1"][rel])
            (work/name).write_bytes(raw)
        # Only original source translation units + new ANTO diagnostic driver
        # with a minimal type alias header. No main.cpp or MCMC compiled.
        compiled=[str(author_git_path/"code"/name) for name in ORIGINAL_CPP]
        exe=base/"phlox_forward"
        command=[compiler,"-std=c++11","-O2","-Wno-deprecated-declarations",
                 "-I",str(shim),"-I",str(author_git_path/"code"),
                 str(DRIVER),*compiled,"-o",str(exe)]
        built=_run_cmd(command,ROOT,timeout=timeout_seconds)
        # 2 replays with the same original author 13-vector parameter START
        first=_run_cmd([str(exe)],work,timeout=timeout_seconds)
        second=_run_cmd([str(exe)],work,timeout=timeout_seconds)
        one=parse_original_driver_output(first.stdout,v42_existing)
        two=parse_original_driver_output(second.stdout,v42_existing)
        # No RNG used in deterministic forward model; ensure replay identical.
        require(one==two,"ORIGINAL_SOURCE_DETERMINISTIC_FORWARD_REPLAY_DRIFT")
        return {
            "schema":"anto_v43_exact_git_pinned_original_PhloxSelection_fixed_start_likelihood",
            "original_paper_DOI":"10.1016/j.cub.2014.07.027",
            "parent_v42_original_full_derived_receipt_replayed_exactly":True,
            "author_original_Git_commit":prev.AUTHOR_COMMIT,
            "original_author_code_unmodified_except_external_type_compatibility_header":True,
            "original_author_missing_nr3_and_Rseeds_not_silently_fabricated":True,
            "source_original_valid_author_input_families_count":371,
            "external_temp_compatibility_nr3_header_sha256":
                hashlib.sha256(COMPAT_HEADER.encode()).hexdigest(),
            "ANTO_authored_single_run_driver_sha256":
                hashlib.sha256(DRIVER.read_bytes()).hexdigest(),
            "source_forward_cpp_translation_units_from_exact_author_commit":
                list(ORIGINAL_CPP),
            "CXX_compiler":compiler,
            "CXX_compiler_version_first_line":
                _run_cmd([compiler,"--version"],ROOT,timeout=10).stdout.splitlines()[0],
            "CXX_compilation_warnings_present":bool(built.stderr.strip()),
            "CXX_compilation_warnings_tail":built.stderr[-1400:],
            "temporary_working_copy_only_all_blank_maternal_padding_lines_removed":
                blank_removed,
            "all_original_nonblank_source_rows_unchanged":True,
            "original_13_start_parameters_from_author_seedSet":[
                .279,.362,.612,.639,.362,0.,.72,.88,.91,.90,3.7,6.3,5.5],
            "model_is_original_2014_forward_one_fixed_set_not_MCMC":True,
            "original_2014_full_75000_step_MCMC_posterior_replayed":False,
            "original_2014_main_paper_full_12024_382_2581_cohort_replayed":False,
            "original_2014_author_EXAMPLE_input_11379_371_2496_replayed":True,
            "source_predicted_probabilities_valid_at_declared_precision":True,
            "independent_CXX_component_log_likelihood_matches_original_calcTotLike":True,
            "two_runs_identical_source_prediction_and_likelihood":True,
            "actual_data_original_forward_output":one,
            "new_empirical_selection_coefficients_inferred":None,
            "natural_ancestral_phenotype_reversion_or_retention_probability_inferred":None,
            "source_original_raw_records_not_redistributed":True,
        }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--author-git-checkout",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--expected",type=Path)
    p.add_argument("--timeout-seconds",type=int,default=140)
    args=p.parse_args()
    result=report(args.author_git_checkout,timeout_seconds=args.timeout_seconds)
    if args.expected and json.loads(args.expected.read_text())!=result:
        raise FixedForwardError("V43_FROZEN_CXX_FIXED_ORIGINAL_SOURCE_LIKELIHOOD_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ANTO_V43_ORIGINAL_FIXED_CPP_SOURCE_REPLAY",json.dumps({
        "original_commit":result["author_original_Git_commit"],
        "expected_parent_source_counts":[11379,371,2496],
        "total":result["actual_data_original_forward_output"]["all_transects"],
        "per_transect":[{k:item[k] for k in (
            "transect","source_model_generations_until_convergence",
            "original_model_log_likelihood",
            "independent_origin_program_flower_component_log_like",
            "independent_origin_program_family_component_log_like")}
            for item in result["actual_data_original_forward_output"]["per_transect"]],
        "fixed_start_not_MCMC":True,
    },sort_keys=True))
    packed=json.dumps(result,sort_keys=True,separators=(",",":"))
    stride=6000
    n=(len(packed)+stride-1)//stride
    for i in range(n):
        print("ANTO_V43_DERIVED_CHUNK",i+1,n,packed[stride*i:stride*(i+1)])


if __name__=="__main__":
    main()
