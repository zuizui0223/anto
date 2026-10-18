"""Tests: source-backed organ contrasts + modular hidden competence theory.

These validations intentionally forbid confusing 2020 C. japonicum organ
measurements with 2026 Ryukyu lineage-specific regain observations.
"""
from dataclasses import replace
from io import BytesIO
from math import log
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import sys
import unittest
import zipfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"scripts"))
from anto_cirsium_organ_evidence_v6 import validate_registry,report
from anto_modular_competence_v6 import (
    Rates,STATES, generator, effective_structural_decay, evolved_distribution,
    observed_flower_distribution, probability_first_reexpression,
    closed_form_retention_then_regulatory_reactivation,synthetic_organ_protection_grid,
)
from probe_cirsium_japonicum_supplement_v6 import inspect_zip,run

EVIDENCE=ROOT/"data/molecular/cirsium_japonicum_organ_evidence_v6.json"


class EmpiricalOrganSourceTest(unittest.TestCase):
    def test_published_contrasts_validated_and_separated(self):
        doc=validate_registry(EVIDENCE)
        self.assertEqual(len(doc["records"]),11)
        self.assertFalse(doc["trait_regain_event_observed"])
        rep=report(EVIDENCE)
        self.assertEqual(rep["flower_vs_leaf_flavonoid_mrna_n"],6)
        self.assertEqual(rep["root_vs_leaf_dfr_mrna_n"],2)
        self.assertFalse(rep["ryukyu_competence_assessed"])
        self.assertFalse(rep["independent_colonization_history_assessed"])

    def test_dfr_denominator_must_not_be_recast_as_flower(self):
        doc=validate_registry(EVIDENCE)
        r=next(x for x in doc["records"] if x["record_id"]=="Park2020_CjDFR_root_isoform_A")
        self.assertEqual((r["numerator_organ"],r["denominator_organ"]),("root","leaf"))
        r["denominator_organ"]="flower"
        with TemporaryDirectory() as temp:
            p=Path(temp)/"fake_source.json"
            p.write_text(json.dumps(doc))
            with self.assertRaises(ValueError):validate_registry(p)

    def test_unmeasured_ryukyu_pigment_not_admitted(self):
        doc=validate_registry(EVIDENCE)
        doc["six_ryukyu_specimens_linked"]=True
        with TemporaryDirectory() as temp:
            p=Path(temp)/"fake_source.json"
            p.write_text(json.dumps(doc))
            with self.assertRaises(ValueError):validate_registry(p)

    def test_contrasts_not_pooled_as_independent_studies(self):
        r=report(EVIDENCE)
        groups={(k["assay"],k["numerator"],k["denominator"]):k["n_distinct_contrasts"]
                for k in r["grouped_by_kind_and_denominator"]}
        self.assertEqual(groups[("qRT-PCR mRNA","flower","leaf")],6)
        self.assertEqual(groups[("qRT-PCR mRNA","root","leaf")],2)
        self.assertEqual(groups[("reported metabolite total","leaf","flower")],1)


class ModularEvolutionModelTest(unittest.TestCase):
    def setUp(self):
        self.r=Rates(.12,0,.1,.05,log(6.0))

    def test_generator_mass_and_positivity(self):
        for u in (0,.25,1):
            Q=generator(self.r,u)
            self.assertTrue(np.allclose(Q.sum(axis=1),0))
            self.assertTrue(np.all(Q-np.diag(np.diag(Q))>=0))
            p=evolved_distribution(self.r,u,"C1_F0",20)
            self.assertAlmostEqual(sum(p.values()),1,places=12)
            self.assertTrue(all(0<=v<=1 for v in p.values()))

    def test_three_different_hidden_white_states_look_identical(self):
        for state in STATES[1:]:
            p={k:float(k==state) for k in STATES}
            observed=observed_flower_distribution(p)
            self.assertEqual(observed,{"pigmented":0,"white":1})

    def test_hidden_functional_states_have_different_futures(self):
        self.assertGreater(probability_first_reexpression(self.r,0,"C1_F0",10),0)
        self.assertEqual(probability_first_reexpression(self.r,0,"C0_F0",10),0)
        self.assertEqual(probability_first_reexpression(self.r,0,"C0_F1",10),0)

    def test_nonfloral_maintenance_can_preserve_coding_but_not_restore_flower_switch(self):
        self.assertAlmostEqual(effective_structural_decay(self.r,1),.02)
        self.assertAlmostEqual(effective_structural_decay(self.r,0),.12)
        p0=closed_form_retention_then_regulatory_reactivation(20,10,1,self.r,0)
        p1=closed_form_retention_then_regulatory_reactivation(20,10,1,self.r,1)
        self.assertGreater(p1,p0)
        zero=replace(self.r,floral_switch_on=0)
        self.assertEqual(closed_form_retention_then_regulatory_reactivation(20,10,1,zero,1),0)

    def test_beta_zero_falsifier(self):
        no_pleiotropy=replace(self.r,nonfloral_constraint_strength=0)
        self.assertAlmostEqual(
            closed_form_retention_then_regulatory_reactivation(20,10,1,no_pleiotropy,0),
            closed_form_retention_then_regulatory_reactivation(20,10,1,no_pleiotropy,1))

    def test_first_hit_handles_alternative_de_novo_path_without_equating_history(self):
        rebuilt=replace(self.r,structural_reconstruction=.09)
        self.assertGreater(probability_first_reexpression(rebuilt,0,"C0_F0",30),0)
        self.assertLessEqual(probability_first_reexpression(rebuilt,0,"C0_F0",30),1)

    def test_illustrative_grid_explicit_not_fitted(self):
        result=synthetic_organ_protection_grid()
        self.assertTrue(result["not_Cirsium_rates"])
        self.assertEqual(len(result["scenario_comparison"]),5)
        self.assertTrue(all(x["nonfloral_use"]>=x["no_nonfloral_use"]
                            for x in result["scenario_comparison"]))

    def test_invalid_rates_rejected(self):
        with self.assertRaises(ValueError): Rates(-.1,0,0,0,0)
        with self.assertRaises(ValueError):generator(self.r,1.5)
        with self.assertRaises(ValueError):probability_first_reexpression(self.r,0,"C1_F1",2)


class SupplementProbeTest(unittest.TestCase):
    def test_offline_never_networks(self):
        self.assertEqual(run()["retrieval_status"],"NOT_REQUESTED")

    def test_small_local_public_fixture_marks_annotations_not_gene_calls(self):
        io=BytesIO()
        with zipfile.ZipFile(io,"w",compression=zipfile.ZIP_DEFLATED) as z:
            z.writestr("ao0c04001_si_002.txt",
                       "geneid\tannotation\nTR001\tdihydroflavonol 4 reductase\nTR002\tcontrol\nTR003\tANS\n")
            workbook=BytesIO()
            with zipfile.ZipFile(workbook,"w",compression=zipfile.ZIP_DEFLATED) as internal:
                internal.writestr("xl/worksheets/sheet1.xml","<worksheet/>")
                internal.writestr("xl/workbook.xml","<workbook/>")
            z.writestr("ao0c04001_si_003.xlsx",workbook.getvalue())
        r=inspect_zip(io.getvalue())
        self.assertEqual(r["annotation_text_scan_status"],"SCANNED_ALL_LINES")
        self.assertEqual(r["potential_gene_annotation_text_hits"],2)
        self.assertEqual(r["focal_ryukyu_evidence"],"NOT_TESTED")
        self.assertEqual(len(r["files_found"]),2)

    def test_decompression_limit(self):
        io=BytesIO()
        with zipfile.ZipFile(io,"w",compression=zipfile.ZIP_DEFLATED) as z:
            z.writestr("ao0c04001_si_002.txt","A"*1024)
        with self.assertRaises(ValueError):inspect_zip(io.getvalue(),max_uncompressed=100)


if __name__=="__main__": unittest.main()
