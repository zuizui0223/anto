"""Adversarial support tests for matched nearest-gap control identifiability."""
import copy
import json
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

from analyze_original_petunieae_phylogeny_v27 import Node
from analyze_petunieae_loo_distance_v31 import tree_distances
import audit_petunieae_gap_identifiability_v34 as m


def star_edges(keys,length):
    return Node(children=[Node(label=k,length=length) for k in keys])


def original_ultrametric_toy():
    # All root->tip lengths equal 5. A and B are each others' closest
    # relatives, and the only size-4 train set with min gap A=10 is C..F.
    root=Node(children=[
        Node(children=[Node(label="A",length=1),Node(label="B",length=1)],length=4),
        Node(children=[Node(label="C",length=1),Node(label="D",length=1)],length=4),
        Node(children=[Node(label="E",length=1),Node(label="F",length=1)],length=4)])
    return root,sorted(root.tips()),["A","B"]


def pairwise_for_labels(ids,override):
    # Positive symmetric default; only distances from focal target matter.
    pairs={a:{b:0.0 if a==b else 30.0 for b in ids} for a in ids}
    for (a,b),val in override.items():
        pairs[a][b]=pairs[b][a]=val
    return pairs


class GapFeasibilityTests(unittest.TestCase):
    def test_frozen_modes_before_source_results(self):
        self.assertEqual(m.CAPS,(8,15))
        self.assertEqual(list(m.TOLERANCES),["abs_1e-5","rel_05pct","rel_10pct","rel_20pct"])
        self.assertEqual(m.TOLERANCES["abs_1e-5"],("absolute",1e-5))
        self.assertEqual(m.TOLERANCES["rel_20pct"],("relative",0.20))
        self.assertEqual(m.EXPECTED_INGROUP,59)

    def test_frozen_tolerance_tuple_json_serialization_is_canonical(self):
        original={"modes":m.TOLERANCES}
        encoded=m.canonical_result(original)
        self.assertEqual(encoded["modes"]["abs_1e-5"],["absolute",1e-5])
        self.assertEqual(encoded,m.canonical_result(encoded))
        self.assertEqual(encoded,json.loads(json.dumps(original)))

    def test_true_complete_ultrametric_clade_has_unique_exact_match(self):
        root,ids,held=original_ultrametric_toy()
        p=tree_distances(root)
        self.assertAlmostEqual(p["A"]["B"],2)
        self.assertAlmostEqual(p["A"]["C"],10)
        x=m.matching_structure(ids,"A",held,p,1e-5)
        self.assertEqual(x["original_training_species_count"],4)
        self.assertEqual(x["original_nearest_training_gap"],10)
        self.assertEqual(x["n_forced_closer_exclusions"],1)
        self.assertEqual(x["n_admissible_training_sets_exact"],1)
        self.assertEqual(x["n_nonidentical_training_sets_exact"],0)
        self.assertEqual(x["n_boundary_candidates"],4)

    def test_gap_band_relaxation_may_or_may_not_create_alternatives(self):
        root,ids,held=original_ultrametric_toy()
        p=tree_distances(root)
        for tol in (1e-5,.5,1.,2.):
            x=m.matching_structure(ids,"A",held,p,tol)
            self.assertEqual(x["n_admissible_training_sets_exact"],1)
        # A third in-clade taxon is already at the boundary in this
        # non-ultrametric example: multiple nonidentical sets are possible.
        ids=list("ABCDEF")
        p=pairwise_for_labels(ids,{
            ("A","B"):8,("A","C"):10,("A","D"):10,
            ("A","E"):11,("A","F"):12,
        })
        x=m.matching_structure(ids,"A",["A","B","C"],p,1e-5)
        self.assertEqual(x["original_training_species_count"],3)
        self.assertEqual(x["n_forced_closer_exclusions"],1)
        self.assertEqual(x["n_boundary_candidates"],2)
        self.assertEqual(x["n_admissible_training_sets_exact"],4)
        self.assertEqual(x["n_nonidentical_training_sets_exact"],3)
        self.assertTrue(x["independent_counterfactual_set_exists"])

    def test_relative_gap_tolerance_in_original_units(self):
        self.assertEqual(m.effective_tolerance(20,("absolute",1e-5)),1e-5)
        self.assertEqual(m.effective_tolerance(20,("relative",.05)),1)
        self.assertEqual(m.effective_tolerance(20,("relative",.1)),2)

    def test_singleton_original_fold_has_no_alternate_training_set(self):
        ids=list("ABCDE")
        p=pairwise_for_labels(ids,{("A","B"):7})
        x=m.matching_structure(ids,"A",["A"],p,1e-5)
        self.assertEqual(x["n_admissible_training_sets_exact"],1)
        self.assertEqual(x["n_nonidentical_training_sets_exact"],0)

    def test_missing_or_duplicate_tips_are_not_imputed(self):
        ids=list("ABCDE")
        p=pairwise_for_labels(ids,{})
        for group in (["A","A"],["B","C"],[]):
            with self.assertRaises(ValueError):
                m.matching_structure(ids,"A",group,p,1e-5)
        with self.assertRaises(ValueError):
            m.matching_structure(["A","A","B","C","D"],"A",["A"],p,1e-5)

    def test_invalid_source_distance_or_geometry_fails(self):
        ids=list("ABCDE")
        p=pairwise_for_labels(ids,{})
        for value in (-1,float("nan"),float("inf")):
            q=copy.deepcopy(p);q["A"]["B"]=value
            with self.assertRaises(ValueError):
                m.matching_structure(ids,"A",["A","B"],q,1e-5)
        q=copy.deepcopy(p);q["A"].pop("C")
        with self.assertRaises(ValueError):
            m.matching_structure(ids,"A",["A","B"],q,1e-5)
        for tol in (-1,float("nan"),float("inf")):
            with self.assertRaises(ValueError):
                m.matching_structure(ids,"A",["A","B"],p,tol)

    def test_eligible_training_sets_monotone_with_band_relaxation(self):
        ids=list("ABCDEF")
        p=pairwise_for_labels(ids,{
            ("A","B"):8,("A","C"):10,("A","D"):10,
            ("A","E"):11,("A","F"):12})
        gaps=[m.matching_structure(ids,"A",["A","B","C"],p,tol)
              ["n_admissible_training_sets_exact"] for tol in (1e-5,.5,1,2,5)]
        self.assertEqual(gaps,sorted(gaps))

    def test_binomial_count_admits_forced_boundary_and_rejects_no_boundary(self):
        self.assertEqual(m.conditional_training_set_count(2,2,3),4)
        self.assertEqual(m.conditional_training_set_count(0,6,3),0)
        self.assertEqual(m.conditional_training_set_count(2,2,5),0)
        for args in ((-1,2,2),(1,-2,2),(1,1,-1)):
            with self.assertRaises(ValueError):
                m.conditional_training_set_count(*args)

    def test_source_fold_aggregation_uses_real_clade_not_independent_tips(self):
        root,ids,held=original_ultrametric_toy()
        p=tree_distances(root)
        folds=[["A","B"],["C","D"],["E","F"]]
        data=m.audit_folds(ids,folds,p)
        self.assertEqual(data["species_scored_once"],6)
        self.assertEqual(data["n_original_clades"],3)
        for mode,cell in data["tolerance_results"].items():
            self.assertEqual(cell["n_targets_with_alternative"],0)
            self.assertEqual(cell["n_original_clades_with_any_alternative"],0)
            self.assertEqual(cell["n_only_original_training_set"],6)
            self.assertEqual(cell["n_not_admitted"],0)
            self.assertEqual(cell["n_targets"],6)

    def test_does_not_drop_nonidentifiable_clades_to_force_positive(self):
        root,ids,held=original_ultrametric_toy()
        p=tree_distances(root)
        report=m.audit_folds(ids,[["A","B"],["C","D"],["E","F"]],p)
        self.assertFalse(report["primary_exact_gap_conditional_comparison_identifiable"])
        self.assertEqual(report["tolerance_results"]["abs_1e-5"]["n_targets_with_alternative"],0)

    def test_any_missing_or_duplicate_fold_target_is_failure(self):
        root,ids,held=original_ultrametric_toy()
        p=tree_distances(root)
        with self.assertRaises(ValueError):
            m.audit_folds(ids,[["A","B"],["A","C"],["D","E","F"]],p)
        with self.assertRaises(ValueError):
            m.audit_folds(ids,[["A","B"],["C","D"]],p)

    def test_biological_observations_cannot_change_structural_geometry(self):
        # The API accepts only the tree+species IDs, not any HPLC/gene data.
        root,ids,held=original_ultrametric_toy()
        p=tree_distances(root)
        before=m.audit_folds(ids,[["A","B"],["C","D"],["E","F"]],p)
        after=m.audit_folds(ids,[["A","B"],["C","D"],["E","F"]],copy.deepcopy(p))
        self.assertEqual(before,after)
        self.assertNotIn("p_value",before)
        self.assertNotIn("gene_function",before)


if __name__=="__main__":
    unittest.main()
