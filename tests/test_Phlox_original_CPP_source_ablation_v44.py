"""ANTO v44 adverse source C++ scenario/relative-likelihood guards (no raw data)."""
from __future__ import annotations

import copy
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_Phlox_original_CPP_source_ablations_v44 as v


class CppPhloxSourcePerturbations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen=v.checked_registry(json.loads(v.REGISTRY.read_text()))
        cls.old=json.loads(v.PARENT.read_text())
        cls.base=cls.old["actual_data_original_forward_output"]

    def test_exact_seven_predeclared_scenarios(self):
        self.assertEqual(v.ORDER,(
            "baseline","fitness_neutral","fitness_uniform_zone0",
            "epistasis_zero","dispersal_half","dispersal_double","assortative_25pct"))
        self.assertEqual(len(self.frozen["scenarios"]),7)
        self.assertTrue(self.frozen["source_model_is_original_2014_forward_not_fitted_MCMC"])

    def test_original_13_parameter_vector_and_true_meaning(self):
        self.assertEqual(len(v.BASELINE),13)
        self.assertEqual(v.BASELINE[5],0)
        self.assertEqual(v.BASELINE[6],.72)
        self.assertEqual(v.BASELINE[7:10],[.88,.91,.9])
        self.assertEqual(v.BASELINE[10:],[3.7,6.3,5.5])

    def test_scenario_mutation_only_declared_positions(self):
        for name,changes in v.CASES.items():
            values=self.frozen["scenarios"][name]["parameter_13_vector"]
            self.assertEqual(
                [i for i,(x,y) in enumerate(zip(v.BASELINE,values)) if x!=y],
                list(changes))
            self.assertEqual(values[:5],v.BASELINE[:5])

    def test_fitness_neutral_retains_baseline_mating_and_dispersal(self):
        x=self.frozen["scenarios"]["fitness_neutral"]["parameter_13_vector"]
        self.assertEqual(x[7:],[1,1,1,1,1,1])
        self.assertEqual(x[5],0)
        self.assertEqual(x[6],.72)

    def test_zone_uniform_keeps_zone_zero_weights_and_source_border(self):
        x=self.frozen["scenarios"]["fitness_uniform_zone0"]["parameter_13_vector"]
        self.assertEqual(x[7:10],x[10:])
        self.assertEqual(x[:7],v.BASELINE[:7])

    def test_epistasis_zero_is_additive_phenotype_model_with_positive_weights(self):
        x=self.frozen["scenarios"]["epistasis_zero"]["parameter_13_vector"]
        self.assertAlmostEqual(1-x[7]+x[8]-x[9],0)
        self.assertAlmostEqual(1-x[10]+x[11]-x[12],0)
        self.assertEqual(x[8],.78)
        self.assertEqual(x[11],8.2)

    def test_geometric_dispersal_not_fraction_of_migrants(self):
        self.assertEqual(self.frozen["scenarios"]["dispersal_half"]["parameter_13_vector"][6],.36)
        self.assertEqual(self.frozen["scenarios"]["dispersal_double"]["parameter_13_vector"][6],1.44)
        self.assertTrue(self.frozen["no_source_biological_fitness_or_allele_frequency_reestimated"])

    def test_assortment_positive_fraction_but_no_new_pollinator_behavior(self):
        x=self.frozen["scenarios"]["assortative_25pct"]["parameter_13_vector"]
        self.assertEqual(x[5],.25)
        self.assertEqual(x[6:],v.BASELINE[6:])

    def test_fake_scenario_registry_data_fails(self):
        s=copy.deepcopy(self.frozen)
        s["scenarios"]["dispersal_double"]["parameter_13_vector"][6]=1.04
        with self.assertRaisesRegex(v.V44Error,"SCENARIO_VALUE_DRIFT"):
            v.checked_registry(s)

    def test_source_v43_driver_reused_with_only_registered_generated_insertions(self):
        text=v.old.DRIVER.read_text()
        result=v.source_derived_cpp(text,self.frozen)
        self.assertIn("V44_REQUIRED_PREDECLARED_SCENARIO",result)
        self.assertIn("V44_UNREGISTERED_SOURCE_SCENARIO",result)
        self.assertIn("V44_SOURCE_WORLD_DID_NOT_REACH_TRUE_CHANGE_TOLERANCE",result)
        self.assertIn("source",str(v.REGISTRY))
        self.assertIn('scenario=="fitness_neutral"',result)
        self.assertIn("const double final_change=world.getMaxChg();",result)
        self.assertIn("calcTotLike",result)
        self.assertNotIn("Mchain",result)
        self.assertNotIn("mcmcStep",result)

    def test_wrong_original_driver_rejected(self):
        wrong=v.old.DRIVER.read_text().replace("std::setprecision(17)","std::setprecision(7)")
        with self.assertRaisesRegex(v.V44Error,"V43_AUTHOR_DIAGNOSTIC_DRIVER_CHANGED"):
            v.source_derived_cpp(wrong,self.frozen)

    def test_five_source_convergence_records(self):
        valid="\n".join(v.TAG_CHANGE+" "+str(i)+" "+str(1e-9)
                        for i in range(1,6))
        self.assertEqual(v.extract_change(valid),{i:1e-9 for i in range(1,6)})

    def test_false_equilibrium_or_missing_deme_rejected(self):
        wrong="\n".join(v.TAG_CHANGE+" "+str(i)+" 0.00000010"
                        for i in range(1,6))
        with self.assertRaisesRegex(v.V44Error,"TRUE_CONVERGENCE"):
            v.extract_change(wrong)
        good="\n".join(v.TAG_CHANGE+" "+str(i)+" 0.0" for i in range(1,5))
        with self.assertRaisesRegex(v.V44Error,"FIVE_TRANSECTS_MISSING"):
            v.extract_change(good)

    def test_exact_original_CXX_baseline_accepted(self):
        c=copy.deepcopy(self.base)
        c["source_true_genotype_change_within_tolerance"]=True
        self.assertTrue(v.verify_base(c,self.base))

    def test_modified_source_likelihood_not_claimed_v43_baseline(self):
        c=copy.deepcopy(self.base)
        c["source_true_genotype_change_within_tolerance"]=True
        c["per_transect"][2]["original_model_log_likelihood"]+=1
        with self.assertRaisesRegex(v.V44Error,"BASELINE_PER_TRANSECT"):
            v.verify_base(c,self.base)

    def test_one_factor_loglike_decomposition_and_sign_counts(self):
        c=copy.deepcopy(self.base)
        for row in c["per_transect"]:
            row["original_model_log_likelihood"]+=2
            row["independent_origin_program_flower_component_log_like"]+=1
            row["independent_origin_program_family_component_log_like"]+=1
        c["all_transects"]["original_total_log_likelihood"]+=10
        c["all_transects"]["original_model_adult_component_log_likelihood"]+=5
        c["all_transects"]["original_model_family_component_log_likelihood"]+=5
        r=v.delta_from_baseline(c,self.base)
        self.assertEqual(r["all_transects_delta_logL"],{"total":10,"adult":5,"offspring":5})
        self.assertEqual(r["n_transects_delta_positive_by_component"],
                         {"total":5,"adult":5,"offspring":5})
        self.assertEqual(r["n_transects_delta_negative_by_component"],
                         {"total":0,"adult":0,"offspring":0})

    def test_component_mismatch_rejected(self):
        c=copy.deepcopy(self.base)
        c["per_transect"][0]["original_model_log_likelihood"]+=4
        c["all_transects"]["original_total_log_likelihood"]+=4
        with self.assertRaisesRegex(v.V44Error,"COMPONENT_DELTA_TOTAL_MISMATCH"):
            v.delta_from_baseline(c,self.base)

    def test_original_main_sample_not_invented(self):
        self.assertTrue(self.frozen["source_32p_11379adult_371mother_2496offspring_unchanged"])
        self.assertEqual(self.base["all_transects"]["n_adults"],11379)
        self.assertEqual(self.base["all_transects"]["n_offspring"],2496)
        self.assertEqual(self.frozen["source_baseline_original_loglike"],-9867.530126894539)

    def test_source_technical_failure_must_remain_not_biological_zero(self):
        with tempfile.TemporaryDirectory() as d:
            fn=Path(d)/"exit3"
            fn.write_text("#!/usr/bin/env sh\nexit 3\n")
            fn.chmod(0o755)
            q=v.source_run_named(fn,Path(d),self.base,2,"fitness_neutral")
            self.assertEqual(q["status"],"BLOCKED_NUMERICAL_OR_ORIGINAL_CODE")
            self.assertEqual(q["failure_category"],"ORIGINAL_CPP_OTHER_NONZERO_EXIT")
            self.assertNotIn("comparison",q)

if __name__=="__main__":
    unittest.main()
