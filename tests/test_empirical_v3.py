"""Prospective negative controls: observable-mechanism identifiability and source gates."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))

from anto_assay_design import emission_matrix, diagnose, run_design_grid
from anto_four_history_design import experiment, comparison, MARKERS
from anto_evidence_gate import load_registry, evaluate, report

REGISTRY=ROOT/'data'/'empirical_sources_v3.json'


class ThreeHistoryAssayTests(unittest.TestCase):
    def test_probabilities(self):
        _, mat=emission_matrix(('foreign_genome','foreign_pigment_haplotype'),missing_probability=.3)
        self.assertTrue(np.allclose(mat.sum(axis=0),1))
        self.assertTrue(np.all((mat>=0)&(mat<=1)))

    def test_colour_only_rank_one(self):
        r=diagnose()
        self.assertEqual(r['rank_known_emissions'],1)
        self.assertAlmostEqual(r['expected_bayes_accuracy'],1/3)
        self.assertAlmostEqual(r['mutual_information_bits'],0)

    def test_one_marker_rank_two(self):
        r=diagnose(('foreign_genome',))
        self.assertEqual(r['rank_known_emissions'],2)
        self.assertFalse(r['mixture_weights_identifiable_if_emissions_known'])

    def test_two_markers_rank_three(self):
        r=diagnose(('foreign_genome','foreign_pigment_haplotype'))
        self.assertEqual(r['rank_known_emissions'],3)
        self.assertGreater(r['expected_bayes_accuracy'],.8)

    def test_uninformative_means_rank_one(self):
        r=diagnose(('foreign_genome','foreign_pigment_haplotype'),.5,.5)
        self.assertEqual(r['rank_known_emissions'],1)
        self.assertAlmostEqual(r['expected_bayes_accuracy'],1/3)

    def test_missingness_reduces_information_not_necessarily_rank(self):
        both=('foreign_genome','foreign_pigment_haplotype')
        full=diagnose(both)
        partly=diagnose(both,missing_probability=.5)
        self.assertLess(partly['expected_bayes_accuracy'],full['expected_bayes_accuracy'])
        self.assertEqual(partly['rank_known_emissions'],3)

    def test_shared_ancestral_haplotype_collapses_distinction(self):
        r=run_design_grid()['shared_ancestral_haplotype']
        self.assertEqual(r['rank_known_emissions'],2)

    def test_invalid_diagnostics_fail(self):
        with self.assertRaises(ValueError):diagnose(('foreign_genome',),1.1,.9)
        with self.assertRaises(ValueError):diagnose(('unknown',))
        with self.assertRaises(ValueError):diagnose(('foreign_genome','foreign_genome'))
        with self.assertRaises(ValueError):diagnose(priors=(0,.5,.5))


class FourHistoryAssayTests(unittest.TestCase):
    def test_colour_only_rank_one(self):
        r=experiment()
        self.assertEqual(r['rank_known_emissions'],1)
        self.assertAlmostEqual(r['bayes_accuracy_known_calibrated_assays'],.25)

    def test_ancestry_alone_cannot_separate_reuse_and_new_mutation(self):
        r=experiment(MARKERS[:2])
        self.assertEqual(r['rank_known_emissions'],3)
        self.assertFalse(r['all_four_mixture_weights_identifiable_if_emissions_known'])

    def test_history_specific_third_measurement_can_increase_rank(self):
        r=experiment(MARKERS)
        self.assertEqual(r['rank_known_emissions'],4)
        self.assertGreater(r['bayes_accuracy_known_calibrated_assays'],.75)

    def test_unvalidated_history_signal_cannot_increase_rank(self):
        r=experiment(MARKERS,collapse_reuse_with_denovo=True)
        self.assertEqual(r['rank_known_emissions'],3)

    def test_weak_assays_are_almost_useless(self):
        a=experiment(MARKERS,.65,.65)
        b=experiment(MARKERS,.9,.95)
        self.assertLess(a['bayes_accuracy_known_calibrated_assays'],
                        b['bayes_accuracy_known_calibrated_assays'])

    def test_invalid_four_way_parameters(self):
        with self.assertRaises(ValueError): experiment(MARKERS,-.1,.9)
        with self.assertRaises(ValueError): experiment((MARKERS[0],MARKERS[0]))
        with self.assertRaises(ValueError): experiment(prior=[.5,.5,0,0])


class EmpiricalGateTests(unittest.TestCase):
    def test_registry_not_yet_jointly_sufficient(self):
        r=report(REGISTRY)
        self.assertEqual(r['integrated_accessibility_systems'],0)
        self.assertEqual(len(r['systems']),4)

    def test_ryukyu_does_not_have_flower_expression(self):
        r=next(r for r in load_registry(REGISTRY) if r['id']=='CIRSIUM_RYUKYU_2026')
        self.assertTrue(r['young_leaf_rna_registered'])
        self.assertFalse(r['flower_stage_rna'])
        self.assertFalse(r['functional_reexpression_test'])
        self.assertEqual(evaluate(r)['admission_tier'],'VISIBLE_LINEAGE_HISTORY_ONLY')

    def test_takaoense_leaf_morphs_are_not_mechanism(self):
        r=next(r for r in load_registry(REGISTRY) if r['id']=='CIRSIUM_TAKAOENSE_2026')
        self.assertTrue(r['same_individual_genomic_flower_link'])
        self.assertFalse(r['floral_pigment_chemistry'])
        self.assertFalse(evaluate(r)['integrated_accessibility_eligible'])

    def test_petunia_is_a_mechanism_not_geographic_rate_control(self):
        r=next(r for r in load_registry(REGISTRY) if r['id']=='PETUNIA_MECHANISTIC_POSITIVE_CONTROL')
        self.assertTrue(evaluate(r)['mechanistic_regain_positive_control'])
        self.assertFalse(evaluate(r)['geographic_history_identified'])

    def test_spurious_functional_flag_is_rejected(self):
        doc=json.loads(REGISTRY.read_text())
        doc['source_records'][0]['functional_reexpression_test']=True
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'invalid.json'
            p.write_text(json.dumps(doc))
            with self.assertRaises(ValueError):load_registry(p)


if __name__=='__main__':
    unittest.main()
