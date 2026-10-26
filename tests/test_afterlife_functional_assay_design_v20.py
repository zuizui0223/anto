"""Strict negative and independent checks for ANTO v20 future assay design."""
import copy
import itertools
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

from anto_functional_assay_identifiability_v20 import panel_from_original_source, sharp_bounds
from anto_functional_assay_bounds_v20 import (
    SOURCE_JOIN, GROUP_SIZE, SOURCE_STATUS, read, function_bound,
    earliest_possible_sign_certificate, build_taxon_ledger, design_receipt
)


def explicit_small_completion(u_total, p_total, u_tested, u_positive,
                              p_tested, p_positive):
    """Enumerate every binary assignment for all remaining small-panel taxa."""
    unobserved_u=u_total-u_tested
    unobserved_p=p_total-p_tested
    outcomes=[]
    for x in itertools.product((0,1),repeat=unobserved_u):
        for y in itertools.product((0,1),repeat=unobserved_p):
            outcomes.append(F(u_positive+sum(x),u_total)
                            - F(p_positive+sum(y),p_total))
    return min(outcomes),max(outcomes)


class FunctionalAssayDesignV20(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.join=read(SOURCE_JOIN)
        cls.result=design_receipt(cls.join)

    def test_original_22_source_taxa_and_source_grain(self):
        r=self.result
        self.assertEqual(r["original_source_taxa"],22)
        self.assertEqual(r["original_source_taxa_by_flower_pigment"],GROUP_SIZE)
        self.assertEqual(r["new_functional_assays_or_regain_events_added"],0)
        self.assertEqual(r["new_multigeneration_reexpression_fitness_experiments"],0)
        self.assertEqual(len(r["taxon_enrollment"]),22)
        self.assertEqual(len({x["taxon"] for x in r["taxon_enrollment"]}),22)
        self.assertTrue(r["entire_clade_phylogenetic_nonindependence_not_resolved"])
        self.assertEqual(r["Cirsium_organismal_empirical_owner"],"zuizui0223/aza3")

    def test_must_never_promote_sequence_annotation_into_function(self):
        rows=self.result["taxon_enrollment"]
        self.assertTrue(all(x["native_flower_joint_function"]=="UNMEASURED" for x in rows))
        self.assertTrue(all(x["floral_induction_response"]=="UNMEASURED" for x in rows))
        self.assertTrue(all(x["transgenerational_stability"]=="UNMEASURED" for x in rows))
        self.assertTrue(all(x["trait_specific_reproduction"]=="UNMEASURED" for x in rows))
        self.assertTrue(all(x["no_new_taxon_function_assay_admitted"] for x in rows))
        self.assertTrue(all("function_present" not in x for x in rows))
        self.assertTrue(self.result["no_cross_angiosperm_rate_or_predictive_clade_holdout_estimated"])

    def test_source_annotation_mask_is_original_three_state_table(self):
        counts=self.result["source_annotation_coverage_not_function"]
        a="DFR_SOURCE_CDS_FEATURE_PRESENT_NOT_FUNCTION"
        b="DFR_SOURCE_RECORD_NO_CDS_FEATURE_NOT_GENE_LOSS"
        c="DFR_NOT_IN_NEW_KT_SOURCE_RANGE_NOT_GENE_LOSS"
        self.assertEqual([counts["UNPIGMENTED"][x] for x in (a,b,c)],[3,2,2])
        self.assertEqual([counts["PIGMENTED"][x] for x in (a,b,c)],[9,0,6])
        self.assertEqual(sum(sum(x.values()) for x in counts.values()),22)
        self.assertEqual(len(SOURCE_STATUS),3)

    def test_no_assays_is_functionally_unidentified_regardless_of_KT_metadata(self):
        current=function_bound(0,0,0,0)
        self.assertEqual(current["unpigmented_minus_pigmented_bound"],["-1","1"])
        self.assertEqual(current["unpigmented_prevalence_bound"],["0","1"])
        self.assertEqual(current["pigmented_prevalence_bound"],["0","1"])
        self.assertEqual(current["finite_panel_sign_decision"],"SIGN_NOT_IDENTIFIED")
        self.assertEqual(self.result["functional_assay_status"],"NO_QUALIFYING_REAL_ASSAYS_ADMITTED")
        self.assertEqual(self.result["identification_if_future_direct_taxon_assays_exist"]["actual_no_new_functional_assays"],current)

    def test_earliest_possible_sign_requires_eight_hypothetical_function_calls(self):
        d=earliest_possible_sign_certificate()
        self.assertEqual(d["minimum_hypothetical_taxon_functional_calls"],8)
        self.assertEqual(d["best_case_allocations_unpigmented_and_pigmented"],[[7,1]])
        self.assertTrue(d["not_a_guaranteed_sample_size_or_power_analysis"])
        self.assertTrue(d["assumes_perfectly_qualified_binary_functional_taxon_assays"])
        # Exhaustive independent integer inequality check: earliest m=8.
        self.assertFalse(any(15*u+7*p>105 for u in range(8)
                             for p in range(16) if u+p<=7 and u<=7))
        self.assertGreater(15*7+7*1,105)

    def test_best_case_positive_direction_strictly_certified(self):
        d=function_bound(7,7,1,0)
        self.assertEqual(d["unpigmented_minus_pigmented_bound"],["1/15","1"])
        self.assertEqual(d["finite_panel_sign_decision"],"UNPIGMENTED_GT_PIGMENTED_CERTIFIED")

    def test_best_case_negative_direction_strictly_certified(self):
        d=function_bound(7,0,1,1)
        self.assertEqual(d["unpigmented_minus_pigmented_bound"],["-1","-1/15"])
        self.assertEqual(d["finite_panel_sign_decision"],"UNPIGMENTED_LT_PIGMENTED_CERTIFIED")

    def test_incomplete_data_cannot_be_force_signed(self):
        d=function_bound(4,2,4,2)
        self.assertEqual(d["unpigmented_minus_pigmented_bound"],["-61/105","61/105"])
        self.assertEqual(d["finite_panel_sign_decision"],"SIGN_NOT_IDENTIFIED")
        self.assertTrue(d["does_not_identify_ancestral_regain_heritability_fitness_or_phylogenetic_causality"])

    def test_full_panel_certifies_exact_descriptive_sample_difference_not_evolution(self):
        d=function_bound(7,5,15,6)
        self.assertEqual(d["unpigmented_minus_pigmented_bound"],["11/35","11/35"])
        self.assertEqual(d["finite_panel_sign_decision"],"UNPIGMENTED_GT_PIGMENTED_CERTIFIED")
        equal=function_bound(7,0,15,0)
        self.assertEqual(equal["unpigmented_minus_pigmented_bound"],["0","0"])
        self.assertEqual(equal["finite_panel_sign_decision"],"EQUALITY_CERTIFIED_FOR_THIS_FINITE_PANEL")

    def test_explicit_binary_completion_oracle_small_finite_populations(self):
        for nu,np in ((2,3),(3,2)):
            for tu in range(nu+1):
                for tp in range(np+1):
                    for pu in range(tu+1):
                        for pp in range(tp+1):
                            d=function_bound(tu,pu,tp,pp,
                              n_unpigmented=nu,n_pigmented=np)
                            a,b=explicit_small_completion(nu,np,tu,pu,tp,pp)
                            self.assertEqual(d["unpigmented_minus_pigmented_bound"],
                                             [str(a),str(b)])

    def test_mathematical_sign_certification_is_not_frequency_or_power_claim(self):
        for key,d in self.result["identification_if_future_direct_taxon_assays_exist"].items():
            self.assertTrue(d["sharp_given_binary_validated_taxon_level_assay"])
            self.assertTrue(d["does_not_identify_ancestral_regain_heritability_fitness_or_phylogenetic_causality"])
            if key.startswith("synthetic_"):
                self.assertNotEqual(d["group_denominators"],[0,0])
        self.assertTrue(self.result["genetic_retention_vs_current_flower_expression_vs_heritable_regain_are_separate_endpoints"])

    def test_direct_taxon_name_relabel_and_function_invention_fail(self):
        d=copy.deepcopy(self.join)
        d["taxon_gene_records"][0]["current_taxon_name"]=d["taxon_gene_records"][1]["current_taxon_name"]
        with self.assertRaisesRegex(ValueError,"distinct"):
            build_taxon_ledger(d)
        d=copy.deepcopy(self.join)
        d["taxon_gene_records"][0]["functional_joint_toolkit"]="FUNCTIONALLY_PROVEN"
        with self.assertRaisesRegex(ValueError,"not validated functional"):
            build_taxon_ledger(d)
        d=copy.deepcopy(self.join)
        d["taxon_gene_records"][0]["direct_flower_reexpression"]="PROVEN"
        with self.assertRaisesRegex(ValueError,"not validated functional"):
            build_taxon_ledger(d)

    def test_original_CDS_feature_and_KT_source_record_must_not_be_swapped(self):
        d=copy.deepcopy(self.join)
        row=next(x for x in d["taxon_gene_records"] if
                 x["DFR_CDS_annotation_in_KT_record"] is False)
        row["DFR_CDS_annotation_in_KT_record"]=None
        with self.assertRaises(ValueError):
            build_taxon_ledger(d)
        d=copy.deepcopy(self.join)
        row=next(x for x in d["taxon_gene_records"] if
                 x["DFR_accession_KT_registered"] is None)
        row["DFR_CDS_annotation_in_KT_record"]=False
        with self.assertRaisesRegex(ValueError,"gene deletion"):
            build_taxon_ledger(d)

    def test_invalid_counts_are_rejected_without_imputation(self):
        bad=((7,8,1,0),(7,7,16,0),(8,7,0,0),(-1,0,0,0),
             (2,3,0,0),(True,0,0,0),(1,0,0,False))
        for args in bad:
            with self.subTest(args=args),self.assertRaises(ValueError):
                function_bound(*args)
        for args in ((0,1),(1,0),(0,0)):
            with self.subTest(args=args),self.assertRaises(ValueError):
                earliest_possible_sign_certificate(*args)

    def test_original_source_row_order_does_not_change_taxon_enrollment(self):
        d=copy.deepcopy(self.join)
        d["taxon_gene_records"].reverse()
        self.assertEqual(build_taxon_ledger(d),self.result["taxon_enrollment"])

    def test_independent_second_v20_model_agrees_on_all_count_combinations(self):
        """Both independently written exact implementations agree for every
        possible measured-positive count combination in the 7/15 panel.

        This is an independent cross-algorithm verification, not additional
        data or an independent real floral regain replication.
        """
        source_panel=panel_from_original_source(self.join)
        our_panel=self.result["taxon_enrollment"]
        self.assertEqual({z["taxon"] for z in source_panel},
                         {z["taxon"] for z in our_panel})
        taxa_u=sorted(x["taxon"] for x in source_panel
                      if x["flower_pigment_category"]=="UNPIGMENTED")
        taxa_p=sorted(x["taxon"] for x in source_panel
                      if x["flower_pigment_category"]=="PIGMENTED")
        self.assertEqual((len(taxa_u),len(taxa_p)),(7,15))
        compared=0
        for tested_u in range(8):
            for positive_u in range(tested_u+1):
                white={name:(i<positive_u)
                       for i,name in enumerate(taxa_u[:tested_u])}
                for tested_p in range(16):
                    for positive_p in range(tested_p+1):
                        labels={**white,**{name:(i<positive_p)
                                            for i,name in enumerate(taxa_p[:tested_p])}}
                        first=function_bound(tested_u,positive_u,tested_p,positive_p)
                        other=sharp_bounds(source_panel,labels,synthetic=True)
                        self.assertEqual(
                            other["difference_unpigmented_minus_pigmented_lower"],
                            F(first["unpigmented_minus_pigmented_bound"][0]))
                        self.assertEqual(
                            other["difference_unpigmented_minus_pigmented_upper"],
                            F(first["unpigmented_minus_pigmented_bound"][1]))
                        compared+=1
        self.assertEqual(compared,4896)

    def test_source_all_registered_accessions_not_an_inherited_assay(self):
        taxa=self.result["taxon_enrollment"]
        self.assertEqual([x["taxon"] for x in taxa],
                         sorted(x["taxon"] for x in taxa))
        ie=next(x for x in taxa if x["taxon"]=="Iochroma edule")
        self.assertEqual(ie["original_flower_pigment_category"],"PIGMENTED")
        self.assertEqual(ie["DFR_new_KT_accession"],"KT898445")
        self.assertEqual(ie["native_ancestral_regain_vs_donor_origin"],"UNRESOLVED")
        no_CDS=[x for x in taxa if x["source_DFR_annotation_class"]==
                "DFR_SOURCE_RECORD_NO_CDS_FEATURE_NOT_GENE_LOSS"]
        self.assertEqual(sorted(x["DFR_new_KT_accession"] for x in no_CDS),
                         ["KT898443","KT898446"])


if __name__=="__main__":
    unittest.main()
