"""Evidence-grain and source-censoring counterfactual tests for ANTO v19.

Unlike a synthetic regain-rate simulator, this validates *real original
source metadata* and protects the distinction between absent deposited
sequence, absent CDS annotation, coding-region claims by original authors,
and actual living floral pathway function.
"""
import json
import sys
import unittest
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from anto_afterlife_empirical_admissibility_v19 import (
    SOURCE_JOIN,REGISTRY,PENSTEMON,IOCHROMA,STATUS,read,analysis,
    contingency,conditional_metadata_reference,record_status,x2,
)


class AfterlifeSourceAdmissibilityV19(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=read(SOURCE_JOIN)
        cls.literature=read(REGISTRY)
        cls.penstemon=read(PENSTEMON)
        cls.iochroma=read(IOCHROMA)
        cls.out=analysis(cls.source,cls.literature,cls.penstemon,cls.iochroma)

    def test_exact_original_22_taxon_three_state_source_table(self):
        self.assertEqual(contingency(self.source),{
            "UNPIGMENTED":dict(zip(STATUS,(3,2,2))),
            "PIGMENTED":dict(zip(STATUS,(9,0,6)))})
        d=self.out["Iochrominae_KT_source_metadata"]
        self.assertEqual(d["original_taxa"],22)
        self.assertEqual(d["source_nucleotide_accessions"],58)
        self.assertEqual(d["DFR_new_KT_accessions"],14)
        self.assertEqual(d["DFR_CDS_features_annotated_in_new_block"],12)
        self.assertEqual(d["DFR_records_with_no_CDS_feature"],2)
        self.assertEqual(d["DFR_no_new_KT_accession_species"],8)

    def test_not_all_missing_cds_are_missing_gene_not_even_one(self):
        r=self.source["taxon_gene_records"]
        nofeature=[x for x in r if x["DFR_CDS_annotation_in_KT_record"] is False]
        not_new=[x for x in r if x["DFR_accession_KT_registered"] is None]
        self.assertEqual({x["DFR_accession_KT_registered"] for x in nofeature},
                         {"KT898443","KT898446"})
        self.assertEqual(len(not_new),8)
        self.assertEqual(len(nofeature),2)
        self.assertTrue(self.out["Iochrominae_KT_source_metadata"][
            "unknown_DFR_new_KT_feature_statuses_do_not_establish_pseudogenization"])
        self.assertTrue(self.out["unmeasured_function_partial_identification"][
            "no_gene_deletion_or_pseudogene_inferred"])

    def test_57_percent_naive_false_loss_calls_conflict_with_original_author(self):
        o=self.out["Iochrominae_KT_source_metadata"]
        self.assertEqual(o["unpigmented_without_new_KT_annotated_DFR"],4)
        self.assertEqual(o["unpigmented_missing_KT_annotated_DFR_fraction"],"4/7")
        self.assertEqual(o[
            "false_negative_if_missing_new_DFR_CDS_feature_treated_as_gene_absence_unpigmented"],4)
        self.assertTrue(o[
            "original_author_reports_intact_three_coding_regions_after_floral_pigment_loss"])
        self.assertTrue(o[
            "author_full_coding_in_original_alignments_not_independently_recomputed"])

    def test_comparing_new_source_annotation_is_not_functional_genetic_contrast(self):
        d=self.out["unmeasured_function_partial_identification"]
        self.assertEqual(self.out["Iochrominae_KT_source_metadata"][
            "observed_annotation_fraction_difference_unpig_minus_pig"],"-6/35")
        self.assertEqual(d["unpigmented_proxy_completion_fraction_range"],["3/7","1"])
        self.assertEqual(d["pigmented_proxy_completion_fraction_range"],["3/5","1"])
        self.assertEqual(d["unpig_minus_pig_proxy_contrast_sharp_bounds"],["-4/7","2/5"])
        self.assertTrue(d["source_proxy_contrast_sign_unidentified"])
        self.assertEqual(d[
            "true_functional_joint_toolkit_fraction_identification_bounds_both_groups"],["0","1"])
        self.assertTrue(d["true_functional_joint_toolkit_unmeasured_in_all_22_by_THIS_KT_SOURCE_AUDIT"])
        self.assertTrue(d["true_group_function_difference_undetermined_without_functional_assay"])

    def test_exact_conditional_metadata_reference_not_biological_significance(self):
        d=self.out["exact_source_metadata_random_label_reference_NOT_evolutionary_p_value"]
        self.assertEqual(d["metadata_observed_pearson_X2_exact"],"33/7")
        self.assertEqual(d["fixed_margins_conditional_tail_exact"],"134/969")
        self.assertEqual(d["enumerated_2_by_3_tables"],21)
        self.assertEqual(d["among_14_deposited_DFR_two_no_CDS_both_unpigmented_exact"],"10/91")
        self.assertTrue(d["conditional_p_is_not_a_phylogenetic_or_evolutionary_test"])
        self.assertGreater(F(d["fixed_margins_conditional_tail_exact"]),F(1,10))

    def test_cross_programme_207_is_not_a_regain_denominator(self):
        d=self.out["cross_programme_evidence_grain"]
        self.assertEqual([x["observations"] for x in d["source_programmes"]],[22,24,161])
        self.assertEqual(d["incompatible_original_observation_units_raw_arithmetic_sum"],207)
        self.assertTrue(d["raw_sum_must_not_be_called_regain_events_or_independent_replicates"])
        self.assertEqual(d["petunia_molecular_mechanism_source_cases"],2)
        self.assertEqual(d["petunia_independent_genus_units"],1)
        self.assertEqual(d["independent_evolutionary_regain_clades_molecularly_source_supported_in_selected_v11_registry"],1)
        self.assertTrue(d["multiple_iochroma_studies_are_within_same_iochrominae_lineage_group"])
        self.assertEqual(d["source_cases_with_trait_specific_qualified_viable_seed_fitness"],0)
        self.assertFalse(d["real_plant_generational_ecological_rate_estimated_by_anto"])
        self.assertFalse(d["observed_species_level_phylogenetic_transition_directly_reconstructed_by_anto"])

    def test_not_a_separate_cirsium_empirical_pipeline(self):
        self.assertEqual(self.out["Cirsium_empirical_owner"],"zuizui0223/aza3")
        self.assertTrue(self.out["no_natural_floral_trait_regain_rate_estimated"])
        self.assertTrue(self.out["no_long_term_fitness_or_adaptation_identified"])
        self.assertTrue(self.out["source_reuse_not_original_author_experiment_or_original_full_alignment_reconstruction"])

    def test_input_duplicate_taxon_rejected(self):
        p=deepcopy(self.source)
        p["taxon_gene_records"][1]["current_taxon_name"]=p["taxon_gene_records"][0]["current_taxon_name"]
        with self.assertRaisesRegex(ValueError,"distinct"):
            contingency(p)

    def test_wrong_source_accession_not_false_gene_loss(self):
        p=deepcopy(self.source)
        p["taxon_gene_records"][0]["CHI_accession"]="KT000000"
        with self.assertRaisesRegex(ValueError,"source|accession"):
            contingency(p)
        p=deepcopy(self.source)
        p["taxon_gene_records"][0]["DFR_accession_KT_registered"]="KT000000"
        with self.assertRaisesRegex(ValueError,"DFR accession"):
            contingency(p)

    def test_taxonomic_rename_is_preserved_not_double_counted(self):
        p=deepcopy(self.source)
        x=next(x for x in p["taxon_gene_records"] if x["taxonomic_join_type"]=="POST2016_COMBINATION_OR_SYNONYM")
        x["taxonomic_join_type"]="DIRECT_BINOMIAL"
        with self.assertRaisesRegex(ValueError,"synonyms"):
            contingency(p)

    def test_status_swaps_fail_even_when_margins_unchanged(self):
        p=deepcopy(self.source)
        a=next(r for r in p["taxon_gene_records"] if r["trait_state_original_fig1"]=="UNPIGMENTED"
               and r["DFR_CDS_annotation_in_KT_record"] is True)
        b=next(r for r in p["taxon_gene_records"] if r["trait_state_original_fig1"]=="PIGMENTED"
               and r["DFR_accession_KT_registered"] is None)
        a["trait_state_original_fig1"]="PIGMENTED"
        b["trait_state_original_fig1"]="UNPIGMENTED"
        with self.assertRaisesRegex(ValueError,"7/15"):
            contingency(p)

    def test_two_metadata_absence_kinds_cannot_be_silently_confused(self):
        a={"DFR_accession_KT_registered":None,
           "DFR_CDS_annotation_in_KT_record":None,
           "DFR_KT_registration_status":"NO_DFR_IN_KT_RANGE_NOT_GENE_ABSENCE"}
        self.assertEqual(record_status(a),STATUS[2])
        a["DFR_CDS_annotation_in_KT_record"]=False
        with self.assertRaisesRegex(ValueError,"gene deletion"):
            record_status(a)
        a["DFR_accession_KT_registered"]="KT898443"
        a["DFR_KT_registration_status"]="RECORDED"
        self.assertEqual(record_status(a),STATUS[1])

    def test_fake_published_mechanistic_independence_rejected(self):
        l=deepcopy(self.literature)
        p=next(x for x in l["cases"] if x["clade"]=="Petunia"
               and x["independent_regain_clade_weight"]==0)
        p["independent_regain_clade_weight"]=1
        with self.assertRaisesRegex(ValueError,"one independent"):
            analysis(self.source,l,self.penstemon,self.iochroma)

    def test_published_donor_cross_not_native_repair(self):
        p=deepcopy(self.penstemon)
        p["reported_source_molecular_results"]["no_direct_natural_pbarbatus_backmutation_reversal_observed"]=False
        with self.assertRaisesRegex(ValueError,"not be natural backmutation"):
            analysis(self.source,self.literature,p,self.iochroma)

    def test_one_backcross_not_independent_reexpression_history(self):
        p=deepcopy(self.iochroma)
        p["not_independent_plant_clade_or_natural_flower_regain"]=False
        with self.assertRaisesRegex(ValueError,"not a natural flower regain"):
            analysis(self.source,self.literature,self.penstemon,p)

    def test_original_f2_count_mutation_fails(self):
        p=deepcopy(self.penstemon)
        p["observed_F2_table"]["red"]["Ff"]+=1
        with self.assertRaisesRegex(ValueError,"matrix drift"):
            analysis(self.source,self.literature,p,self.iochroma)

    def test_original_2016_coding_alignment_never_reidentified_from_KT_CDS(self):
        p=deepcopy(self.source)
        p["functional_biochemistry_or_expression_verified_here"]=True
        with self.assertRaisesRegex(ValueError,"claim boundaries"):
            contingency(p)
        p=deepcopy(self.source)
        p["taxon_gene_records"][0]["functional_joint_toolkit"]="FUNCTIONAL"
        with self.assertRaisesRegex(ValueError,"not validated functional"):
            contingency(p)


if __name__=="__main__":
    unittest.main()
