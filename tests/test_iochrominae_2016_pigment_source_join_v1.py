"""Figure-1 flower colour × original NCBI accession crosswalk tests."""
import sys
from pathlib import Path
from copy import deepcopy
import unittest
import json

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_iochrominae_2016_pigment_source_join_v1 import calculate,load,GENES,TRAITS,ANNOTATIONS

class IochrominaeTraitSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g=load(GENES);cls.t=load(TRAITS);cls.a=load(ANNOTATIONS)
    def test_frozen_original_join_exact(self):
        expected=load(ROOT/"results/iochrominae_2016_trait_KT_source_join_20261019.json")
        self.assertEqual(calculate(self.g,self.t,self.a),expected)
        self.assertEqual(expected["n_species"],22)
        c=expected["original_source_counts"]
        self.assertEqual(c["UNPIGMENTED"]["species"],7)
        self.assertEqual(c["PIGMENTED"]["species"],15)
        self.assertEqual(c["UNPIGMENTED"]["CHI_KT_records"],7)
        self.assertEqual(c["UNPIGMENTED"]["F3H_KT_records"],7)
        self.assertEqual(c["UNPIGMENTED"]["DFR_KT_records"],5)
        self.assertEqual(c["PIGMENTED"]["DFR_KT_records"],9)
    def test_annotated_cds_is_not_recoverability_or_function(self):
        d=calculate(self.g,self.t,self.a)
        white=d["original_source_counts"]["UNPIGMENTED"]
        colored=d["original_source_counts"]["PIGMENTED"]
        self.assertEqual((white["DFR_KT_records_without_CDS_feature"],
                          colored["DFR_KT_records_without_CDS_feature"]),(2,0))
        self.assertTrue(d["CDS_feature_missing_not_enzyme_loss"])
        self.assertFalse(d["functional_biochemistry_or_expression_verified_here"])
        self.assertFalse(d["independent_cross_clade_regain_rate_estimated"])
    def test_Iochroma_edule_retained_records_is_not_a_new_causal_reversal(self):
        x=next(x for x in calculate(self.g,self.t,self.a)["taxon_gene_records"]
               if x["current_taxon_name"]=="Iochroma edule")
        self.assertEqual((x["CHI_accession"],x["F3H_accession"],x["DFR_accession_KT_registered"]),
                         ("KT898423","KT898401","KT898445"))
        self.assertEqual(x["historical_regain_evidence"],"PUBLISHED_PHYLOGENETIC_INFERENCE")
        self.assertEqual(x["functional_joint_toolkit"],"UNMEASURED")
    def test_wrong_taxon_concept_refused(self):
        x=deepcopy(self.t);x["source_taxa"][0]["current_genbank_taxon"]="Iochroma madeup"
        with self.assertRaisesRegex(ValueError,"crosswalk changed"):calculate(self.g,x,self.a)
    def test_pigment_reclassification_does_not_pass(self):
        x=deepcopy(self.t);x["source_taxa"][0]["floral_anthocyanin_status_from_2016_figure"]="UNPIGMENTED" if x["source_taxa"][0]["floral_anthocyanin_status_from_2016_figure"]=="PIGMENTED" else "PIGMENTED"
        with self.assertRaisesRegex(ValueError,"biological status counts"):calculate(self.g,x,self.a)
    def test_missing_KT_record_does_not_become_gene_deletion(self):
        x=next(x for x in calculate(self.g,self.t,self.a)["taxon_gene_records"]
               if x["current_taxon_name"]=="Iochroma loxense")
        self.assertIsNone(x["DFR_accession_KT_registered"])
        self.assertEqual(x["DFR_KT_registration_status"],"NO_DFR_IN_KT_RANGE_NOT_GENE_ABSENCE")
        self.assertIsNone(x["DFR_CDS_annotation_in_KT_record"])
    def test_exact_source_record_CDS_identity_guard(self):
        x=deepcopy(self.a);x["source_records_without_CDS_annotation"][0]["accession"]="KT000000"
        with self.assertRaisesRegex(ValueError,"CDS-annotation"):calculate(self.g,self.t,x)
if __name__=="__main__":unittest.main()
