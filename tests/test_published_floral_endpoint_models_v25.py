"""Source-published evidence tests for ANTO v25 endpoint-specific prediction.

No original 2019 individual qPCR/HPLC was retrieved and no 2016 gene
function, individual tracking or actual ancestral trait regain is
experimentally established by this published-table reanalysis.
"""
import json
import sys
import unittest
from copy import deepcopy
from decimal import Decimal
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_published_floral_endpoint_models_v25 import (
    SOURCE,RAW_GATE,V24_TABLE,
    EXPECTED_MODELS,EXPECTED_CORRELATIONS,EXPECTED_DIMENSIONS,
    source_verification,conditional_akaike,original_tissue_correlation,
    report,load
)


class OriginalPublishedFloralEndpointModelsV25(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=load(SOURCE)
        cls.gate=load(RAW_GATE)
        cls.v24=load(V24_TABLE)
        cls.output=report(cls.source,cls.gate,cls.v24)

    def test_original_two_paper_provenance_and_original_retrieval_abstention(self):
        self.assertTrue(source_verification(self.source,self.gate,self.v24))
        self.assertEqual(self.source["source_2016"]["doi"],
                         "10.1016/j.phytochem.2016.05.007")
        self.assertEqual(self.source["source_2019"]["doi"],"10.1002/dvdy.82")
        self.assertFalse(self.gate["original_RAR_bytes_verified"])
        self.assertEqual(self.gate["raw_gene_expression_sample_rows_parsed"],0)
        self.assertTrue(self.output["no_original_2019_individual_HPLC_qPCR_data_replayed"])
        self.assertTrue(self.output["not_independent_plant_clades_or_same_original_individuals"])

    def test_two_original_published_five_best_ou_pgls_model_tables_are_exact(self):
        d=self.source["source_2019"]["reported_models"]
        self.assertEqual([x["response"] for x in d],["anthocyanidins","flavonols"])
        for response in d:
            got=[(x["model"],x["aic"],x["delta_aic"],x["npar"])
                 for x in response["models"]]
            self.assertEqual(got,EXPECTED_MODELS[response["response"]])
            self.assertEqual(len(got),5)
            self.assertEqual(got[0][0],"Chi")

    def test_aic_chi_model_best_for_both_quantitative_outputs_only(self):
        output=self.output["published_quantitative_outcomes"]
        for name in ("anthocyanidins","flavonols"):
            x=output[name]
            self.assertEqual(x["published_best_model"],"Chi")
            self.assertEqual(x["published_runner_up"],"F3h")
            self.assertEqual(x["n_listed_candidates"],5)
            self.assertTrue(x["full_candidate_set_completeness_unverified"])
            self.assertTrue(x["model_fits_not_independently_recomputed"])
            self.assertAlmostEqual(sum(z[
                "within_five_candidate_conditional_Akaike_weight"]
                for z in x["all_five"]),1,places=10)

    def test_two_original_printed_AIC_gaps_and_relative_evidence_ratios(self):
        a=self.output["published_quantitative_outcomes"]["anthocyanidins"]
        b=self.output["published_quantitative_outcomes"]["flavonols"]
        self.assertEqual(a["best_vs_runner_up_delta_AIC"],"12.2")
        self.assertEqual(b["best_vs_runner_up_delta_AIC"],"9.9")
        self.assertAlmostEqual(a["relative_best_vs_runner_up_AIC_evidence_ratio"],
                               445.857770,places=5)
        self.assertAlmostEqual(b["relative_best_vs_runner_up_AIC_evidence_ratio"],
                               141.174964,places=5)
        self.assertAlmostEqual(a["conditional_Chi_weight_among_only_five_printed"],
                               0.997758915569,places=10)
        self.assertAlmostEqual(b["conditional_Chi_weight_among_only_five_printed"],
                               0.992963340356,places=10)
        for x in (a,b):
            lo,hi=x["relative_best_vs_runner_up_ratio_interval_for_rounded_AIC_plusminus_0p05"]
            self.assertLess(lo,x["relative_best_vs_runner_up_AIC_evidence_ratio"])
            self.assertGreater(hi,x["relative_best_vs_runner_up_AIC_evidence_ratio"])

    def test_2016_real_original_9x9_matrix_symmetric_signed_not_rounded_to_binary(self):
        d=self.source["source_2016"]
        self.assertEqual(d["original_individuals_in_published_phylogenetic_correlation"],26)
        self.assertEqual(d["original_species_in_published_phylogenetic_correlation"],9)
        self.assertEqual(d["reported_phylogenetic_lambda"],"0.342")
        self.assertEqual(tuple(d["dimension_labels"]),EXPECTED_DIMENSIONS)
        self.assertEqual(tuple(map(tuple,d["original_printed_correlation_matrix"])),
                         EXPECTED_CORRELATIONS)
        self.assertEqual(len(d["original_printed_correlation_matrix"]),9)
        self.assertEqual(len(d["original_printed_correlation_matrix"][0]),9)

    def test_2016_original_petal_branch_correlations_and_leaf_tissue_specificity(self):
        r=self.output["published_tissue_metabolite_correlations"]
        p=r["selected_original_printed_r"]
        self.assertEqual(p["within_petal_di_antho_vs_di_flav"],"17/25")
        self.assertEqual(p["within_petal_tri_antho_vs_tri_flav"],"13/25")
        self.assertEqual(p["within_petal_mono_flav_vs_di_flav"],"39/50")
        self.assertEqual(p["within_leaf_mono_flav_vs_di_flav"],"24/25")
        self.assertEqual(p["cross_petal_di_antho_vs_leaf_di_flav"],"3/5")
        self.assertEqual(p["cross_petal_di_antho_vs_leaf_mono_flav"],"27/50")
        self.assertEqual(r["unique_within_tissue_pairs"],18)
        self.assertEqual(r["unique_between_tissue_pairs"],18)
        self.assertEqual(r["within_tissue_pairs_abs_r_at_least_0p5"],4)
        self.assertEqual(r["between_tissue_pairs_abs_r_at_least_0p5"],2)
        self.assertTrue(r["one_published_matrix_not_36_independent_plant_experiments"])
        self.assertTrue(r["flower_to_leaf_cross_tissue_r_not_zero_or_universally_strong"])

    def test_aic_support_is_within_five_listed_not_all_model_posterior(self):
        for name,d in self.output["published_quantitative_outcomes"].items():
            self.assertEqual(len(d["all_five"]),5)
            self.assertTrue(d["full_candidate_set_completeness_unverified"])
            self.assertTrue(all(z[
                "model_weight_is_not_gene_causality_or_heldout_prediction"]
                for z in d["all_five"]))
            self.assertNotIn("probability_ChI_causally_regulates_trait",d)

    def test_original_colour_classification_cannot_substitute_for_abundance_AIC(self):
        d=self.output["source_reported_white_fixed_vs_white_morph_gene_expression_class_contrast"]
        self.assertFalse(d["new_held_out_classification_model_fitted_by_ANTO"])
        self.assertTrue(d["does_not_equal_quantitative_pigment_PGLS_ranking"])
        self.assertIn("DFR/ANS",d["original_paper_claim"])
        self.assertEqual(self.output["horizon_target_separation"][
            "multi_generation_trait_reexpression"],"NOT OBSERVED OR ESTIMATED")

    def test_alter_2019_model_entry_aic_or_delta_or_identity_is_rejected(self):
        for field,new in (
            ("aic","77.4"),("delta_aic","12.3"),
            ("npar",20),("model","Dfr"),("genes",["DFR"])):
            p=deepcopy(self.source)
            p["source_2019"]["reported_models"][0]["models"][1][field]=new
            with self.subTest(field=field),self.assertRaises(ValueError):
                source_verification(p,self.gate,self.v24)

    def test_alter_2016_correlation_coefficient_or_observation_unit_rejected(self):
        for target,value in (
           ("corr","0.69"),("n",27),("species",10),("lambda","0.75"),
           ("gene_labels","leaf_mono_anthocyanidin")):
            d=deepcopy(self.source)
            early=d["source_2016"]
            if target=="corr":
                early["original_printed_correlation_matrix"][3][4]=value
                early["original_printed_correlation_matrix"][4][3]=value
            elif target=="n":early["original_individuals_in_published_phylogenetic_correlation"]=value
            elif target=="species":early["original_species_in_published_phylogenetic_correlation"]=value
            elif target=="lambda":early["reported_phylogenetic_lambda"]=value
            else:early["dimension_labels"][2]=value
            with self.subTest(target=target),self.assertRaises(ValueError):
                source_verification(d,self.gate,self.v24)

    def test_no_raw_data_can_be_promoted_from_published_author_model(self):
        d=deepcopy(self.gate)
        d["original_RAR_bytes_verified"]=True
        with self.assertRaisesRegex(ValueError,"raw-material availability"):
            source_verification(self.source,d,self.v24)
        d=deepcopy(self.source)
        d["scope_control"]["no_original_2019_hplc_qpcr_raw_reads"]=False
        with self.assertRaisesRegex(ValueError,"Biological data status"):
            source_verification(d,self.gate,self.v24)
        d=deepcopy(self.v24)
        d["individually_measured_2019_gene_expression_or_HPLC_obtained"]=True
        with self.assertRaisesRegex(ValueError,"raw-material availability"):
            source_verification(self.source,self.gate,d)

    def test_published_2016_tissues_not_the_2019_exact_same_plant_samples(self):
        r=self.output
        self.assertTrue(r["not_independent_plant_clades_or_same_original_individuals"])
        self.assertEqual(r["published_tissue_metabolite_correlations"][
            "published_correlation_species"],9)
        self.assertTrue(r["source_reported_white_fixed_vs_white_morph_gene_expression_class_contrast"][
            "does_not_equal_quantitative_pigment_PGLS_ranking"])
        self.assertTrue(r["no_natural_colour_return_or_ancestral_toolkit_function_confirmed"])
        self.assertTrue(r["not_posterior_gene_causal_probability_or_cross_clade_effect"])

    def test_attempt_to_claim_flower_function_or_ancestral_regain_rejected(self):
        d=deepcopy(self.source)
        for key in ("no_new_direct_flower_joint_functional_assays",
                    "no_natural_ancestral_colour_regain_rate",
                    "no_between_source_causal_transport"):
            p=deepcopy(d)
            p["scope_control"][key]=False
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,"Biological data status"):
                source_verification(p,self.gate,self.v24)

    def test_original_table_paper_mislabeled_rejected(self):
        d=deepcopy(self.source)
        d["source_2019"]["doi"]="10.1093/molbev/msy117"
        with self.assertRaisesRegex(ValueError,"Original Table2"):
            source_verification(d,self.gate,self.v24)
        d=deepcopy(self.source)
        d["source_2016"]["doi"]="10.1002/dvdy.82"
        with self.assertRaisesRegex(ValueError,"Original Table2"):
            source_verification(d,self.gate,self.v24)

    def test_full_source_receipt_and_optional_exact_freeze(self):
        r=self.output
        self.assertEqual(r["schema"],
                         "anto_v25_original_literature_endpoint_separation_flavonoid_evidence")
        self.assertTrue(r["no_original_2019_individual_HPLC_qPCR_data_replayed"])
        self.assertEqual(r["Cirsium_organismal_empirical_owner"],"zuizui0223/aza3")
        json.dumps(r,sort_keys=True)
        receipt=ROOT/"results/iochrominae_published_pathway_endpoint_v25_20261028.json"
        if receipt.is_file():
            self.assertEqual(r,load(receipt))


if __name__=="__main__":
    unittest.main()
