"""Adversarial checks for v43 original fixed-start C++ source diagnostic."""
import copy
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import run_original_PhloxSelection_fixed_forward_v43 as v


def fake_valid_log():
    rows=[]
    for t,(p,ad,m,k) in v.EXPECTED_WORLD_COUNTS.items():
        la=-float(ad);lf=-float(k);ll=la+lf
        probs="0.25 0.25 0.25 0.25"
        rows.append(
            f"{v.TAG_WORLD} {t} {p} {ad} {m} {k} 234 {ll} {la} {lf} 0.0 {probs} {probs}"
        )
    rows.append(f"{v.TAG_TOTAL} 32 11379 371 2496 -13875 -11379 -2496 0")
    return "\n".join(rows)


class OriginalFixedCPPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent=json.loads(v.V42_RESULT.read_text())

    def test_eight_pinned_original_forward_translation_units_not_original_MCMC(self):
        self.assertEqual(len(v.ORIGINAL_CPP),8)
        self.assertNotIn("mcmc.cpp",v.ORIGINAL_CPP)
        self.assertNotIn("main.cpp",v.ORIGINAL_CPP)
        self.assertIn("likelihood.cpp",v.ORIGINAL_CPP)
        self.assertIn("World.cpp",v.ORIGINAL_CPP)

    def test_source_specific_compat_nr3_only_pure_aliases(self):
        x=v.COMPAT_HEADER
        self.assertIn("typedef std::vector<double> VecDoub;",x)
        self.assertIn("typedef const VecDoub& VecDoub_I;",x)
        self.assertIn("typedef double Doub;",x)
        self.assertNotIn("NRvector",x)
        self.assertNotIn("Mchain",x)

    def test_original_source_13_parameters_frozen_and_unfitted(self):
        self.assertEqual(v.EXPECTED_SOURCE_START,
                         [.279,.362,.612,.639,.362,0,.72,.88,.91,.9,3.7,6.3,5.5])
        self.assertEqual(self.parent["mother_source_family_row_total"],371)
        self.assertEqual(self.parent["source_offspring_phenotype_total"],2496)

    def test_strip_blank_family_padding_only(self):
        original=(b"population\tcolor\tDB\tDR\tLB\tLR\n"
                  b"315\t3\t1\t0\t2\t0\n"
                  b"315\t3\t0\t3\t0\t1\n"
                  b"\t\t\t\t\t\n\n")
        out,n=v.strip_only_blank_original_family_padding(original)
        self.assertEqual(n,2)
        self.assertEqual(out,(b"population\tcolor\tDB\tDR\tLB\tLR\n"
                              b"315\t3\t1\t0\t2\t0\n"
                              b"315\t3\t0\t3\t0\t1\n"))

    def test_nonblank_bad_family_row_is_not_dropped_to_manufacture_success(self):
        original=(b"population\tcolor\tDB\tDR\tLB\tLR\n"
                  b"315\tbad\t1\t0\t2\t0\n"
                  b"\t\t\t\t\t\n")
        out,n=v.strip_only_blank_original_family_padding(original)
        self.assertEqual(n,1)
        self.assertIn(b"bad",out)
        with self.assertRaises(v.prev.OriginalSourceError):
            v.prev.parse_fam(out,{315})

    def test_bad_source_family_header_must_fail(self):
        with self.assertRaisesRegex(v.FixedForwardError,"AUTHOR_FAMILY_HEADER"):
            v.strip_only_blank_original_family_padding(
                b"mother\tcolor\tDB\tDR\tLB\tLR\n315\t1\t1\t0\t0\t0\n"
            )

    def test_five_transects_original_fixed_source_pop_and_family_counts(self):
        self.assertEqual(v.EXPECTED_WORLD_COUNTS,{
            1:(9,4562,121,649),2:(7,1744,87,630),3:(6,2300,56,407),
            4:(6,1424,78,586),5:(4,1349,29,224)})
        self.assertEqual(sum(x[1] for x in v.EXPECTED_WORLD_COUNTS.values()),11379)
        self.assertEqual(sum(x[2] for x in v.EXPECTED_WORLD_COUNTS.values()),371)
        self.assertEqual(sum(x[3] for x in v.EXPECTED_WORLD_COUNTS.values()),2496)

    def test_parse_valid_five_source_components_and_total(self):
        r=v.parse_original_driver_output(fake_valid_log(),self.parent)
        self.assertEqual(len(r["per_transect"]),5)
        self.assertEqual(r["all_transects"]["n_adults"],11379)
        self.assertEqual(r["all_transects"]["n_mothers"],371)
        self.assertEqual(r["all_transects"]["n_offspring"],2496)
        self.assertEqual(r["all_transects"]["original_total_log_likelihood"],-13875.)

    def test_missing_transect_or_duplicate_source_result_fails(self):
        good=fake_valid_log()
        for mutated in (
            "\n".join(x for x in good.splitlines() if "SOURCE_TRANSECT 5 " not in x),
            good+"\n"+next(x for x in good.splitlines() if "SOURCE_TRANSECT 2 " in x),
        ):
            with self.assertRaises(v.FixedForwardError):
                v.parse_original_driver_output(mutated,self.parent)

    def test_likelihood_part_component_drift_rejected(self):
        good=fake_valid_log()
        wrong=good.replace(
            "ANTO_V43_ORIGINAL_SOURCE_TOTAL 32 11379 371 2496 -13875",
            "ANTO_V43_ORIGINAL_SOURCE_TOTAL 32 11379 371 2496 -13876")
        with self.assertRaisesRegex(v.FixedForwardError,"LIKELIHOOD"):
            v.parse_original_driver_output(wrong,self.parent)

    def test_prediction_row_sum_nonunit_and_nan_fail(self):
        good=fake_valid_log()
        for bad in ("0.95 0.25 0.25 0.25","nan 0.25 0.25 0.25"):
            wrong=good.replace("0.25 0.25 0.25 0.25",bad,1)
            with self.assertRaises(v.FixedForwardError):
                v.parse_original_driver_output(wrong,self.parent)

    def test_author_source_likelihood_log_not_positive(self):
        good=fake_valid_log()
        original=next(x for x in good.splitlines() if x.startswith(
            v.TAG_WORLD+" 1 "))
        fields=original.split()
        self.assertEqual(float(fields[7]),-5211.0)
        fields[7]="1.0"
        wrong=good.replace(original," ".join(fields))
        with self.assertRaises(v.FixedForwardError):
            v.parse_original_driver_output(wrong,self.parent)

    def test_no_fake_mother_family_padded_row_count(self):
        good=fake_valid_log()
        wrong=good.replace(" 1424 78 586 "," 1424 80 586 ",1)
        with self.assertRaisesRegex(v.FixedForwardError,"DENOMINATOR"):
            v.parse_original_driver_output(wrong,self.parent)

    def test_no_false_2014_main_12k_adult_as_fixed_author_example(self):
        good=fake_valid_log()
        wrong=good.replace("SOURCE_TOTAL 32 11379 371 2496",
                           "SOURCE_TOTAL 32 12024 382 2581")
        with self.assertRaisesRegex(v.FixedForwardError,"DENOMINATOR"):
            v.parse_original_driver_output(wrong,self.parent)

    def test_model_no_start_tuning_after_outcome(self):
        driver=v.DRIVER.read_text()
        self.assertIn("ORIGINAL_SEEDSET_VALUES_DRIFT",driver)
        self.assertIn("calcTotLike",driver)
        self.assertIn("MakeLikelihood",driver)
        self.assertIn("SOURCE_ORIGINAL_LIKELIHOOD_COMPONENTS_DISAGREE",driver)
        self.assertNotIn("Mchain",driver)
        self.assertNotIn("mcmcStep",driver)

    def test_replay_parent_source_digest_never_omitted(self):
        src=v.V42_RESULT
        self.assertTrue(src.exists())
        self.assertTrue(self.parent["original_author_source"][
            "all_source_digests_verified_from_pinned_original_checkout"])
        self.assertEqual(self.parent["source_populations_with_no_example_maternal_records"],
                         [330,332])

    def test_technical_compiler_timeout_is_not_biological_zero(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(v.FixedForwardError):
                v._run_cmd(["python3","-c","import time;time.sleep(4)"],
                           Path(d),timeout=.01)

if __name__=="__main__":
    unittest.main()
