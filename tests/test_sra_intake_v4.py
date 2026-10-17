"""Regression checks on real public ENA run metadata and aza3 morph crosswalk."""
import csv
import json
from collections import Counter
from copy import deepcopy
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_sra_intake import load_rows, validate_takaoense, validate_platform, audit

MANIFEST=ROOT/"data/public_sra/prjna1311153_run_manifest_20261017.csv"
MORPH=ROOT/"data/public_sra/takaoense_aza3_morph_crosswalk_20261017.csv"
PLATFORM=ROOT/"data/public_sra/prjna1311153_ena_platform_fingerprints_20261017.json"


def save_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class SourceSpecimenCrosswalkTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=load_rows(MANIFEST)

    def test_25_distinct_run_biosample_experiments(self):
        self.assertEqual(len(self.rows),25)
        for key in ("run_accession","experiment_accession","biosample_accession","specimen_code_from_experiment_title"):
            self.assertEqual(len(set(r[key] for r in self.rows)),25)

    def test_ryukyu_six_exactly_not_colour_observations(self):
        subset=[r for r in self.rows if r["taxon_as_experiment_titled"] in
                ("Cirsium brevicaule","Cirsium irumtiense")]
        self.assertEqual(len(subset),6)
        self.assertEqual(Counter(r["taxon_as_experiment_titled"] for r in subset),
                         Counter({"Cirsium brevicaule":3,"Cirsium irumtiense":3}))
        self.assertEqual({r["specimen_code_from_experiment_title"] for r in subset},
                         {"AMM-4163","YRN-4166","OKN-4295","MMY-4078","ISG-4111","YGN-4296"})
        self.assertTrue(all(r["individual_flower_colour"]=="NOT_MEASURED_IN_THIS_DATASET" for r in subset))

    def test_aza3_six_morphs_join_by_run_and_sample_and_specimen(self):
        crosswalk=validate_takaoense(MORPH,self.rows)
        self.assertEqual(len(crosswalk),6)
        self.assertEqual(Counter(r["morph_code_from_aza3"] for r in crosswalk),
                         Counter({"W":3,"BP":3}))

    def test_declared_md5s_and_novaseq_miniseq_mismatch_preserved(self):
        platform=validate_platform(PLATFORM,self.rows)
        self.assertEqual(platform["ena_instrument_models"]["Illumina MiniSeq"],25)
        self.assertIn("NovaSeq",platform["paper_methods_instrument"])
        self.assertTrue(platform["metadata_discrepancy"])

    def test_offline_report_does_not_claim_regain(self):
        r=audit(MANIFEST,MORPH,PLATFORM)
        self.assertEqual(r["ryukyu_n_runs"],6)
        self.assertEqual(r["takaoense_verified_morph_crosswalk_n"],6)
        self.assertFalse(r["floral_reexpression_demonstrated"])
        self.assertFalse(r["ryukyu_individual_pigment_assayed"])

    def test_forged_white_pigment_from_taxon_is_rejected(self):
        rows=deepcopy(self.rows)
        rows[0]["individual_flower_colour"]="white"
        with TemporaryDirectory() as d:
            p=Path(d)/"fake.csv"
            save_csv(p,rows)
            with self.assertRaises(ValueError):load_rows(p)

    def test_swapped_biosample_is_rejected_by_morph_join(self):
        rows=deepcopy(self.rows)
        r=next(x for x in rows if x["run_accession"]=="SRR35152717")
        wrong=next(x for x in rows if x["run_accession"]=="SRR35152718")
        r["biosample_accession"], wrong["biosample_accession"]=wrong["biosample_accession"],r["biosample_accession"]
        with self.assertRaises(ValueError): validate_takaoense(MORPH,rows)

    def test_wrong_specimen_title_must_fail(self):
        rows=deepcopy(self.rows)
        rows[0]["specimen_code_from_experiment_title"]="WRONG-999"
        with TemporaryDirectory() as d:
            p=Path(d)/"fake.csv"
            save_csv(p,rows)
            with self.assertRaises(ValueError):load_rows(p)

    def test_unsupported_plant_tissue_assumption_rejected(self):
        rows=deepcopy(self.rows)
        rows[0]["tissue_from_project_description"]="flower corolla"
        with TemporaryDirectory() as d:
            p=Path(d)/"fake.csv"
            save_csv(p,rows)
            with self.assertRaises(ValueError):load_rows(p)

    def test_changed_run_count_rejected(self):
        with TemporaryDirectory() as d:
            p=Path(d)/"fake.csv"
            save_csv(p,self.rows[:-1])
            with self.assertRaises(ValueError):load_rows(p)


if __name__=="__main__":
    unittest.main()
