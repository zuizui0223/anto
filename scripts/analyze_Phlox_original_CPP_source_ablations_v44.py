"""ANTO v44: evaluate pre-frozen 7 one-factor changes in original C++ Phlox model.

Reuse Git-SHA-pinned original 2014 PhloxSelection C++ (8 unmodified forward
translation units), the ANTO-authored v43 diagnostic driver, and original
source EXAMPLE records. Compare actual author 5-transect log likelihoods for
adult four-color counts vs mother-conditioned child four-color counts. An
unoptimized start-parameter sensitivity is NOT a newly estimated fitness or
migration parameter; no MCMC, P-values or independent natural replication.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import subprocess
import tempfile
from pathlib import Path

import run_original_PhloxSelection_fixed_forward_v43 as old
ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v44_source_cpp_fixed_scenario_vectors_20261117.json"
PARENT=ROOT/"results/anto_v43_original_Phlox_fixed_forward_likelihood_20261117.json"
BASELINE=[.279,.362,.612,.639,.362,0,.72,.88,.91,.90,3.7,6.3,5.5]
CASES={
 "baseline":{},
 "fitness_neutral":{7:1,8:1,9:1,10:1,11:1,12:1},
 "fitness_uniform_zone0":{10:.88,11:.91,12:.90},
 "epistasis_zero":{8:.78,11:8.2},
 "dispersal_half":{6:.36},
 "dispersal_double":{6:1.44},
 "assortative_25pct":{5:.25},
}
ORDER=tuple(CASES)
TAG_CHANGE="ANTO_V44_ACTUAL_ORIGINAL_WORLD_GENOTYPE_CHANGE"


class V44Error(ValueError):
    pass


def require(ok,why):
    if not ok: raise V44Error(why)


def checked_registry(source):
    require(source.get("schema")==
        "anto_v44_original_2014_author_cpp_predeclared_parameter_ablation_vectors",
        "V44_REGISTERED_SCENARIO_SOURCE_SCHEMA_DRIFT")
    require(source.get("source_original_git_sha")==old.prev.AUTHOR_COMMIT
            and source.get("source_original_repo")==old.prev.AUTHOR_REPO
            and source.get("parent_v43_source_start_parameter_vector")==BASELINE
            and source.get("scenario_ids_fixed_in_order")==list(ORDER)
            and source.get("source_baseline_original_loglike")==-9867.530126894539
            and source.get("model_demography_convergence_original_tolerance")==1e-8,
            "V44_SOURCE_OR_PRIOR_AUTHOR_START_NOT_FIXED")
    require(source.get("model_original_parameter_indices")=={
        "border_1_to_5":[0,1,2,3,4],
        "assortative_mating_a":5,
        "source_geometric_pollen_dispersal_m":6,
        "zone0_fitness_DB_DR_LR":[7,8,9],
        "zone1_fitness_DB_DR_LR":[10,11,12],
    },"ORIGINAL_CPP_PARAMETER_INDEX_IDENTITY_DRIFT")
    require(set(source.get("scenarios",{}))==set(ORDER),
            "SOURCE_SCENARIO_SET_DRIFT")
    for name,changes in CASES.items():
        r=source["scenarios"][name]
        expected=[changes.get(i,x) for i,x in enumerate(BASELINE)]
        require(r["parameter_13_vector"]==expected and
                r["changed_zero_based_parameter_positions"]==list(changes)
                and len(expected)==13 and all(math.isfinite(x) for x in expected),
                "SOURCE_PRE_FROZEN_SCENARIO_VALUE_DRIFT_"+name)
        require(all(x>0 for x in expected[7:]) and
                0<=expected[5]<=1 and expected[6]>0,
                "SOURCE_SCENARIO_PARAMETER_OUTSIDE_BIOLOGICALLY_MEANINGFUL_DOMAIN")
    require(source["source_model_is_original_2014_forward_not_fitted_MCMC"] is True
            and source["source_32p_11379adult_371mother_2496offspring_unchanged"] is True
            and source["no_source_biological_fitness_or_allele_frequency_reestimated"] is True,
            "V44_SCIENTIFIC_SCOPE_GUARD_DRIFT")
    return source


def source_derived_cpp(base,source):
    require(hashlib.sha256(base.encode()).hexdigest()==
            json.loads(PARENT.read_text())["ANTO_authored_single_run_driver_sha256"],
            "V43_AUTHOR_DIAGNOSTIC_DRIVER_CHANGED")
    require(base.count("int main() {")==1 and
            base.count("std::vector<double> start=original_start();")==1 and
            base.count("int generations=(int)world.output_Gens();")==1,
            "V44_ALLOWED_DRIVER_INSERTION_POINTS_CHANGED")
    require(base.count("if (generations==params.maxGens)")==1,
            "V44_AUTHOR_ORIGINAL_FORWARD_CONVERGENCE_GUARD_DRIFT")
    source=checked_registry(source)
    cases=[]
    for name in ORDER:
        values=source["scenarios"][name]["parameter_13_vector"]
        seq=", ".join(format(float(x),".17g") for x in values)
        cases.append(
            ('if' if not cases else 'else if')+
            ' (scenario=="'+name+'") start = {'+seq+'};'
        )
    cpp=base.replace("int main() {","int main(int argc, char** argv) {",1)
    injection=(
        'std::vector<double> start=original_start();\n'
        '    if (argc!=2) throw std::runtime_error("V44_REQUIRED_PREDECLARED_SCENARIO");\n'
        '    const std::string scenario(argv[1]);\n'
        '    '+"\n    ".join(cases)+'\n'
        '    else throw std::runtime_error("V44_UNREGISTERED_SOURCE_SCENARIO");'
    )
    cpp=cpp.replace("std::vector<double> start=original_start();",injection,1)
    needle="int generations=(int)world.output_Gens();"
    tag=(
        needle+"\n"
        '        const double final_change=world.getMaxChg();\n'
        '        if (!std::isfinite(final_change) || final_change>params.maxChange+1e-14)\n'
        '            throw std::runtime_error("V44_SOURCE_WORLD_DID_NOT_REACH_TRUE_CHANGE_TOLERANCE");\n'
        '        std::cout<<"'+TAG_CHANGE+' "<<params.wName[i]<<" "<<final_change<<std::endl;'
    )
    cpp=cpp.replace(needle,tag,1)
    require("mcmcStep" not in cpp and "Mchain" not in cpp,
            "V44_UNADMITTED_MCMC_DRIVER_MODIFICATION")
    return cpp


def extract_change(stdout):
    changes={}
    for line in stdout.splitlines():
        if not line.startswith(TAG_CHANGE+" "):
            continue
        cells=line[len(TAG_CHANGE)+1:].split()
        require(len(cells)==2,"V44_GENOTYPE_CONVERGENCE_RECORD_MALFORMED")
        try:
            transect=int(cells[0]); value=float(cells[1])
        except ValueError as exc:
            raise V44Error("V44_BAD_CONVERGENCE_NUMERIC") from exc
        require(transect in range(1,6) and transect not in changes
                and math.isfinite(value) and 0<=value<=1e-8+1e-14,
                "V44_SOURCE_TRUE_CONVERGENCE_NOT_CONFIRMED")
        changes[transect]=value
    require(set(changes)==set(range(1,6)),
            "V44_SOURCE_TRUE_CONVERGENCE_FIVE_TRANSECTS_MISSING")
    return changes


def source_run(exe,work,old_receipt,timeout_seconds):
    try:
        run=subprocess.run([str(exe)],cwd=work,text=True,
                           stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                           timeout=timeout_seconds,check=False)
    except subprocess.TimeoutExpired:
        return {"status":"BLOCKED_NUMERICAL_OR_ORIGINAL_CODE",
                "failure_category":"ORIGINAL_CPP_NUMERICAL_OR_RUNTIME_TIMEOUT"}
    if run.returncode:
        text=run.stderr+" "+run.stdout[-350:]
        if "V44_SOURCE_WORLD_DID_NOT_REACH_TRUE_CHANGE_TOLERANCE" in text:
            reason="ORIGINAL_CPP_STABILITY_STOP_WITH_NONZERO_FINAL_CHANGE"
        elif "FIXED_SOURCE_WORLD_HIT_MAX_GENERATIONS" in text:
            reason="ORIGINAL_CPP_MAX_GENERATIONS_REACHED"
        elif "MODEL_NOT_NORMALIZED" in text or "INVALID_PROBABILITY" in text:
            reason="ORIGINAL_CPP_PROBABILITY_INVALID"
        elif "ORIGINAL_SOURCE_NONFINITE" in text:
            reason="ORIGINAL_CPP_NONFINITE_LOG_LIKELIHOOD"
        elif run.returncode<0:
            reason="ORIGINAL_CPP_PROCESS_SIGNAL"
        else:
            reason="ORIGINAL_CPP_OTHER_NONZERO_EXIT"
        return {"status":"BLOCKED_NUMERICAL_OR_ORIGINAL_CODE",
                "failure_category":reason,
                "process_exit_code":run.returncode}
    try:
        value=old.parse_original_driver_output(run.stdout,old_receipt)
        changes=extract_change(run.stdout)
    except (old.FixedForwardError,V44Error) as exc:
        return {"status":"BLOCKED_NUMERICAL_OR_ORIGINAL_CODE",
                "failure_category":"SOURCE_OR_NUMERICAL_OUTPUT_GATE_"+type(exc).__name__,
                "gate_message":str(exc)}
    for x in value["per_transect"]:
        x["final_max_genotype_frequency_change"]=changes[x["transect"]]
    value["source_true_genotype_change_within_tolerance"]=True
    return {"status":"ADMITTED_ORIGINAL_FORWARD_SOURCE", "output":value}


def delta_from_baseline(value,base):
    rows=[]
    for v,b in zip(value["per_transect"],base["per_transect"]):
        require(v["transect"]==b["transect"],
                "V44_TRANSECT_ALIGNMENT_MISMATCH")
        rows.append({
            "transect":v["transect"],
            "Delta_logL_total_vs_original_START":
                v["original_model_log_likelihood"]-b["original_model_log_likelihood"],
            "Delta_logL_adult_vs_original_START":
                v["independent_origin_program_flower_component_log_like"]-
                b["independent_origin_program_flower_component_log_like"],
            "Delta_logL_offspring_vs_original_START":
                v["independent_origin_program_family_component_log_like"]-
                b["independent_origin_program_family_component_log_like"],
        })
    o=value["all_transects"];p=base["all_transects"]
    sums={}
    for child,parent,name in (
        ("original_total_log_likelihood","original_total_log_likelihood","total"),
        ("original_model_adult_component_log_likelihood",
         "original_model_adult_component_log_likelihood","adult"),
        ("original_model_family_component_log_likelihood",
         "original_model_family_component_log_likelihood","offspring")
    ):
        sums[name]=o[child]-p[parent]
    require(math.isclose(sums["total"],sums["adult"]+sums["offspring"],
                         rel_tol=0,abs_tol=1e-6),
            "V44_COMPONENT_DELTA_TOTAL_MISMATCH")
    for name,label in (("total","total"),("adult","adult"),("offspring","offspring")):
        key="Delta_logL_"+label+"_vs_original_START"
        require(math.isclose(sums[name],sum(v[key] for v in rows),
                             rel_tol=0,abs_tol=1e-6),
                "V44_PER_TRANSECT_DELTA_SUM_MISMATCH")
    positive={kind:sum(row["Delta_logL_"+kind+"_vs_original_START"]>1e-7
                       for row in rows) for kind in sums}
    negative={kind:sum(row["Delta_logL_"+kind+"_vs_original_START"]< -1e-7
                       for row in rows) for kind in sums}
    return {"all_transects_delta_logL":sums,
            "per_transect_delta_logL":rows,
            "n_transects_delta_positive_by_component":positive,
            "n_transects_delta_negative_by_component":negative,
            "component_deltas_exactly_add_to_total_within_numeric_tolerance":True}


def verify_base(value,old_data):
    require(value["source_true_genotype_change_within_tolerance"],
            "V43_BASELINE_FALSE_EQUILIBRIUM")
    now=value["all_transects"]; frozen=old_data["all_transects"]
    for key in ("original_total_log_likelihood",
                "original_model_adult_component_log_likelihood",
                "original_model_family_component_log_likelihood"):
        require(math.isclose(now[key],frozen[key],rel_tol=0,abs_tol=1e-7),
                "V43_FIXED_BASELINE_TOTAL_LIKELIHOOD_DRIFT")
    for a,b in zip(value["per_transect"],old_data["per_transect"]):
        require(a["transect"]==b["transect"] and
                a["source_model_generations_until_convergence"]==
                    b["source_model_generations_until_convergence"],
                "V43_FIXED_BASELINE_SOURCE_GENERATIONS_DRIFT")
        for key in ("original_model_log_likelihood",
                    "independent_origin_program_flower_component_log_like",
                    "independent_origin_program_family_component_log_like"):
            require(math.isclose(a[key],b[key],rel_tol=0,abs_tol=1e-7),
                    "V43_FIXED_BASELINE_PER_TRANSECT_LIKELIHOOD_DRIFT")
    return True


def report(source_dir,timeout_seconds=160):
    registry=checked_registry(json.loads(REGISTRY.read_text()))
    v43=json.loads(PARENT.read_text())
    manifest=json.loads(old.prev.MANIFEST.read_text())
    v41=json.loads(old.prev.V41.read_text())
    v42=json.loads(old.V42_RESULT.read_text())
    current=old.prev.report(source_dir,manifest,v41)
    require(current==v42,
            "V44_PARENT_V42_VERIFIED_AUTHOR_INPUT_DRIFT")
    require(v43["author_original_Git_commit"]==old.prev.AUTHOR_COMMIT
            and v43["two_runs_identical_source_prediction_and_likelihood"] is True,
            "V44_PARENT_V43_FROZEN_RESULT_DRIFT")
    require(not (source_dir/"code/nr3.h").exists()
            and not (source_dir/"Rseeds.txt").exists(),
            "V44_AUTHOR_ORIGINAL_DEPENDENCY_SET_CHANGED")
    source_driver=old.DRIVER.read_text()
    injected=source_derived_cpp(source_driver,registry)
    with tempfile.TemporaryDirectory(prefix="anto_v44_") as tmp:
        tmp=Path(tmp)
        work=tmp/"work";shim=tmp/"shim"
        work.mkdir();shim.mkdir()
        (shim/"nr3.h").write_text(old.COMPAT_HEADER)
        blank={}
        for t in range(1,6):
            for suffix in ("Phen","Fam"):
                fname=f"inPhlox_{suffix}_T{t}.data"
                rel="inputFiles/"+fname
                raw=old.prev.digest_file(source_dir/rel,manifest[
                     "all_author_files_fixed_git_blob_sha1"][rel])
                if suffix=="Fam":
                    cleaned,n_blank=old.strip_only_blank_original_family_padding(raw)
                    blank[fname]=n_blank
                    # Source author reader cannot skip padded blank mother rows
                    original_families=old.prev.parse_fam(cleaned,{
                        x["population"] for x in old.prev.parse_phen(old.prev.digest_file(
                            source_dir/("inputFiles/inPhlox_Phen_T"+str(t)+".data"),
                            manifest["all_author_files_fixed_git_blob_sha1"][
                                "inputFiles/inPhlox_Phen_T"+str(t)+".data"]))})
                    require(len(original_families)==old.EXPECTED_WORLD_COUNTS[t][2],
                            "V44_MOTHER_ONLY_BLANK_PADDING_SOURCE_ROW_MISMATCH")
                    raw=cleaned
                (work/fname).write_bytes(raw)
            for name in (f"init_Phlox_IF{t}.txt",):
                rel="inputFiles/"+name
                (work/name).write_bytes(old.prev.digest_file(source_dir/rel,manifest[
                     "all_author_files_fixed_git_blob_sha1"][rel]))
        for name in ("inPhlox_nuissance.params","seedSet.params"):
            rel="parameterFiles/"+name
            (work/name).write_bytes(old.prev.digest_file(source_dir/rel,manifest[
                "all_author_files_fixed_git_blob_sha1"][rel]))
        require(blank=={
            "inPhlox_Fam_T1.data":0,"inPhlox_Fam_T2.data":0,
            "inPhlox_Fam_T3.data":0,"inPhlox_Fam_T4.data":299,
            "inPhlox_Fam_T5.data":0},
            "V44_AUTHOR_SOURCE_PADDING_NOT_THE_2014_EXAMPLE")
        generated=tmp/"v44_fixed.cpp"
        generated.write_text(injected)
        exe=tmp/"forward_scenarios"
        compile_args=["g++","-std=c++11","-O2","-Wno-deprecated-declarations",
                      "-I",str(shim),"-I",str(source_dir/"code"),
                      str(generated),*[
                          str(source_dir/"code"/name) for name in old.ORIGINAL_CPP
                      ],"-o",str(exe)]
        build=old._run_cmd(compile_args,ROOT,timeout=timeout_seconds)
        cases={}
        for case in ORDER:
            def run_one():
                # Call separate executable process, report technical failure
                return source_run(tmp/"forward_scenarios",work,v42,timeout_seconds)
            # passing the only permitted case name is handled below through
            # a wrapper in the "run_one" closure, never by source output
            first=source_run_named(exe,work,v42,timeout_seconds,case)
            if first["status"]!="ADMITTED_ORIGINAL_FORWARD_SOURCE":
                if case=="baseline":
                    raise V44Error("ORIGINAL_BASELINE_FIXED_SOURCE_BLOCKED_"+first["failure_category"])
                cases[case]={**first,"parameters":registry["scenarios"][case],
                             "numerical_failure_not_biological_zero":True}
                continue
            second=source_run_named(exe,work,v42,timeout_seconds,case)
            require(second==first,
                    "V44_SCENARIO_FORWARD_TWO_IDENTICAL_SOURCE_RUNS_DISAGREE_"+case)
            v=first["output"]
            if case=="baseline":
                verify_base(v,v43["actual_data_original_forward_output"])
            cases[case]={"status":"ADMITTED_ORIGINAL_FORWARD_SOURCE",
                         "parameters":registry["scenarios"][case],
                         "output":v,"two_forward_runs_identical":True,
                         "source_true_change_tolerance_verified":True}
        base=cases["baseline"]["output"]
        for case in ORDER:
            if cases[case]["status"]=="ADMITTED_ORIGINAL_FORWARD_SOURCE":
                cases[case]["comparison"]=delta_from_baseline(cases[case]["output"],base)
                if case=="baseline":
                    require(math.isclose(cases[case]["comparison"]["all_transects_delta_logL"]["total"],
                                         0,rel_tol=0,abs_tol=1e-10),
                            "BASELINE_NOT_ZERO_RELATIVE_TO_ITSELF")
        admitted=[c for c in ORDER if cases[c]["status"]=="ADMITTED_ORIGINAL_FORWARD_SOURCE"]
        blocked=[c for c in ORDER if cases[c]["status"]!="ADMITTED_ORIGINAL_FORWARD_SOURCE"]
        return {
            "schema":"anto_v44_original_cpp_fixed_start_seven_source_parameter_ablation",
            "analysis_type":"EXPLORATORY_FIXED_START_SOURCE_CPP_SENSITIVITY_NOT_MODEL_SELECTION",
            "author_original_commit":old.prev.AUTHOR_COMMIT,
            "author_original_source_20_blobs_and_v42_receipt_replayed":True,
            "parent_v43_original_cpp_fixed_baseline_likelihood_reproduced":True,
            "parent_v43_actual_total_logL":v43["actual_data_original_forward_output"][
                "all_transects"]["original_total_log_likelihood"],
            "source_original_cpp_eight_forward_units_unchanged":list(old.ORIGINAL_CPP),
            "source_original_CXX_compatibility_shim_sha256":
                hashlib.sha256(old.COMPAT_HEADER.encode()).hexdigest(),
            "source_ANTO_v43_driver_sha256":
                hashlib.sha256(source_driver.encode()).hexdigest(),
            "temporary_scenario_driver_sha256":
                hashlib.sha256(injected.encode()).hexdigest(),
            "original_source_32p_11379adults_371families_2496offspring":True,
            "original_2014_published_full_12024_382_2581_used":False,
            "fixed_scenario_order":list(ORDER),
            "source_parameter_vectors_not_fitted":[registry["scenarios"][name][
                "parameter_13_vector"] for name in ORDER],
            "n_source_scenarios_admitted":len(admitted),
            "n_source_scenarios_technically_blocked":len(blocked),
            "scenario_names_admitted":admitted,
            "scenario_names_blocked":blocked,
            "source_2014_original_CXX_compile_warnings_present":bool(build.stderr.strip()),
            "source_2014_original_CXX_compile_warnings_tail":build.stderr[-1200:],
            "scenarios":cases,
            "n_independent_original_phylogenetic_studies_added":0,
            "actual_empirical_selection_or_migration_estimated":False,
            "original_2014_75000_step_MCMC_performed":False,
            "new_ancestral_flower_colour_regain_or_allele_frequency_estimated":False,
            "conditional_likelihood_changes_are_not_causal_mechanism_significance":True,
            "no_original_source_raw_rows_redistributed":True,
        }


def source_run_named(exe,work,old_receipt,timeout_seconds,name):
    try:
        p=subprocess.run([str(exe),name],cwd=work,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                         text=True,timeout=timeout_seconds,check=False)
    except subprocess.TimeoutExpired:
        return {"status":"BLOCKED_NUMERICAL_OR_ORIGINAL_CODE",
                "failure_category":"ORIGINAL_CPP_NUMERICAL_OR_RUNTIME_TIMEOUT"}
    if p.returncode:
        text=p.stderr+" "+p.stdout[-400:]
        if "V44_SOURCE_WORLD_DID_NOT_REACH_TRUE_CHANGE_TOLERANCE" in text:
            reason="ORIGINAL_CPP_STABILITY_STOP_WITH_NONZERO_FINAL_CHANGE"
        elif "FIXED_SOURCE_WORLD_HIT_MAX_GENERATIONS" in text:
            reason="ORIGINAL_CPP_MAX_GENERATIONS_REACHED"
        elif "MODEL_NOT_NORMALIZED" in text or "INVALID_PROBABILITY" in text:
            reason="ORIGINAL_CPP_PROBABILITY_INVALID"
        elif "ORIGINAL_SOURCE_NONFINITE" in text:
            reason="ORIGINAL_CPP_NONFINITE_LOG_LIKELIHOOD"
        elif p.returncode<0:
            reason="ORIGINAL_CPP_PROCESS_SIGNAL"
        else:
            reason="ORIGINAL_CPP_OTHER_NONZERO_EXIT"
        return {"status":"BLOCKED_NUMERICAL_OR_ORIGINAL_CODE",
                "failure_category":reason,"process_exit_code":p.returncode}
    try:
        result=old.parse_original_driver_output(p.stdout,old_receipt)
        changes=extract_change(p.stdout)
    except (old.FixedForwardError,V44Error) as exc:
        return {"status":"BLOCKED_NUMERICAL_OR_ORIGINAL_CODE",
                "failure_category":"ORIGINAL_SOURCE_OUTPUT_GATES_"+type(exc).__name__,
                "gate_message":str(exc)}
    for x in result["per_transect"]:
        x["final_max_genotype_frequency_change"]=changes[x["transect"]]
    result["source_true_genotype_change_within_tolerance"]=True
    return {"status":"ADMITTED_ORIGINAL_FORWARD_SOURCE","output":result}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--author-git-checkout",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--expected",type=Path)
    p.add_argument("--timeout-seconds",type=int,default=160)
    args=p.parse_args()
    result=report(args.author_git_checkout,timeout_seconds=args.timeout_seconds)
    if args.expected and json.loads(args.expected.read_text())!=result:
        raise V44Error("V44_FROZEN_SOURCE_ABLATION_RESULTS_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    compact={
       name:({
          "status":x["status"],
          "difference":x["comparison"]["all_transects_delta_logL"],
          "n_transects_positive":x["comparison"]["n_transects_delta_positive_by_component"],
       } if x["status"]=="ADMITTED_ORIGINAL_FORWARD_SOURCE" else {
          "status":x["status"],"reason":x["failure_category"]}
       ) for name,x in result["scenarios"].items()}
    print("ANTO_V44_SOURCE_2014_FORWARD_ABLATION_SUMMARY",json.dumps(
          {"scenarios":compact,"source_numerical_baseline_reproduced":True},sort_keys=True))
    packed=json.dumps(result,sort_keys=True,separators=(",",":"))
    step=6000
    n=(len(packed)+step-1)//step
    for i in range(n):
        print("ANTO_V44_DERIVED_CHUNK",i+1,n,packed[i*step:(i+1)*step])


if __name__=="__main__":
    main()
