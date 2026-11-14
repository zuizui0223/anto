"""Adversarial v30 nearest patristic-path and predictor novelty/heldout risk checks."""
from __future__ import annotations

import copy
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_petunieae_patristic_prediction_v30 as m
from analyze_original_petunieae_phylogeny_v27 import parse_newick, TreeGateError
from analyze_original_petunieae_2023_species_prediction_v26 import (
    MODELS, ridge_predict
)


def synthetic_tree(n=31):
    codes=[f"XX{i:02d}" for i in range(n)]
    leaves=[f"{v}:{1+(i%5)/10}" for i,v in enumerate(codes)]
    tree="(BROW:5,("+",".join(leaves)+"):2);"
    return parse_newick(tree),codes


def synthetic_data(codes):
    names=sorted(set(g for xs in MODELS.values() for g in xs))
    out={}
    for i,taxon in enumerate(codes):
        gene={g:math.log1p(1+(i%7)*.2+j*.3) for j,g in enumerate(names)}
        out[taxon]={
            "genes":gene,
            "responses":{
                "anthocyanidins_log1p_total":math.log1p(1+i%11)+.14*gene["FLS"],
                "flavonols_log1p_total":math.log1p(1+i%13)+.13*gene["FLS"],
            },
        }
    return out


def synthetic_records(cap=5):
    tree,codes=synthetic_tree()
    data=synthetic_data(codes)
    pair=m.patristic_pairs(tree)
    folds=[codes[i:i+cap] for i in range(0,31,cap)]
    records=m.source_distances_in_folds(data,folds,pair)
    return tree,codes,data,folds,records


class V30SourceAdverse(unittest.TestCase):
    def test_protocol_all_caps_and_genes_frozen(self):
        self.assertEqual(m.FINE_AND_COARSE_31,(2,3,4,5,6,8,15))
        self.assertEqual(m.COARSE_59,(8,15))
        self.assertEqual(m.BRANCH_GENES,("F3primeH","F35H","FLS"))
        self.assertEqual(m.PAT_UNIT,
            "published_2023_Newick_branch_length_units_not_calibrated_years")
        self.assertEqual(len(MODELS),6)
        self.assertEqual(m.RIDGE_LAMBDA,10)

    def test_exact_shortest_known_newick_path_and_root_stem_ignored(self):
        text="((AAAA:1,BBBB:2):3,"+",".join(f"X{i:03d}:0.5" for i in range(18))+"):900;"
        pair=m.patristic_pairs(parse_newick(text))
        self.assertEqual(len(pair),20)
        self.assertAlmostEqual(pair["AAAA"]["BBBB"],3)
        self.assertAlmostEqual(pair["AAAA"]["X000"],4.5)
        self.assertAlmostEqual(pair["BBBB"]["X000"],5.5)
        self.assertEqual(pair["AAAA"]["AAAA"],0)
        self.assertEqual(pair["X001"]["AAAA"],pair["AAAA"]["X001"])
        self.assertAlmostEqual(pair["X001"]["X002"],1)

    def test_missing_edges_do_not_get_imputed_to_zero(self):
        bad="((AAAA,BBBB:2):3,"+",".join(f"X{i:03d}:0.5" for i in range(18))+");"
        with self.assertRaisesRegex(TreeGateError,"SOURCE_BRANCH_LENGTH_NOT_AVAILABLE"):
            m.patristic_pairs(parse_newick(bad))
        zero="("+",".join(f"X{i:03d}:0" for i in range(20))+");"
        with self.assertRaisesRegex(TreeGateError,"PATRISTIC_ALL_ZERO_SOURCE_TREE"):
            m.patristic_pairs(parse_newick(zero))

    def test_bad_tree_tip_duplicates_and_negative_lengths_rejected(self):
        dup="("+",".join(["AAAA:1","AAAA:2"]+
                        [f"X{i:03d}:1" for i in range(18)])+");"
        with self.assertRaisesRegex(TreeGateError,"NEWICK_EMPTY_OR_DUPLICATE_TIPS"):
            parse_newick(dup)
        neg="("+",".join([f"X{i:03d}:"+("-1" if i==2 else "1")
                            for i in range(20)])+");"
        with self.assertRaisesRegex(TreeGateError,"NEWICK_NEGATIVE_OR_NONFINITE_LENGTH"):
            parse_newick(neg)

    def test_training_novelty_and_test_outcome_never_used(self):
        _,codes,data,folds,records=synthetic_records()
        train=sorted(set(codes)-set(folds[0]))
        target=folds[0][0]
        dist=m.training_only_nearest_expression(target,data,train)
        self.assertGreaterEqual(dist,0)
        changed=copy.deepcopy(data)
        changed[target]["responses"]["anthocyanidins_log1p_total"]=99999
        self.assertEqual(m.training_only_nearest_expression(target,data,train),
                         m.training_only_nearest_expression(target,changed,train))
        with self.assertRaisesRegex(TreeGateError,"MOLECULAR_MODULE_CHANGED"):
            m.training_only_nearest_expression(target,data,train,genes=("CHI",))
        with self.assertRaisesRegex(TreeGateError,"MOLECULAR_MODULE_CHANGED"):
            m.training_only_nearest_expression(target,data,train+[target])

    def test_exact_species_coverage_and_observation_only_prediction(self):
        tree,codes,data,folds,records=synthetic_records()
        self.assertEqual(len(records),31)
        self.assertEqual(set(t["species"] for t in records),set(codes))
        for r in records:
            self.assertEqual(set(r["responses"]),set(m.ENDPOINTS))
            self.assertTrue(math.isfinite(r["patristic"]))
            self.assertTrue(math.isfinite(r["expression_novelty"]))
            self.assertGreaterEqual(r["patristic"],0)
        wrong=folds+[folds[0]]
        with self.assertRaisesRegex(TreeGateError,"PREDICTION_DENOMINATOR_DRIFT"):
            m.source_distances_in_folds(data,wrong,m.patristic_pairs(tree))

    def test_fixed_distance_only_tiers_and_no_target_selection(self):
        _,codes,data,folds,records=synthetic_records()
        tiers,_=m.tier_assignments(records)
        self.assertEqual(sorted(len(z) for z in tiers),[10,10,11])
        self.assertEqual(len(set(x["species"] for t in tiers for x in t)),31)
        control=copy.deepcopy(records)
        for row in control:
            for e in m.ENDPOINTS:
                row["responses"][e]["residuals"][m.BRANCH_MODEL]=-999999
                row["responses"][e]["sq_branch"]=999999
        new,_=m.tier_assignments(control)
        self.assertEqual([[r["species"] for r in t] for t in tiers],
                         [[r["species"] for r in t] for t in new])

    def test_59_original_source_tiers(self):
        reps=[{"species":f"X{i:03d}","patristic":float(i)}
              for i in range(59)]
        tiers,_=m.tier_assignments(reps)
        self.assertEqual(sorted(len(z) for z in tiers),[19,20,20])
        self.assertEqual(sum(len(x) for x in tiers),59)

    def test_midrank_and_spearman_tie_handling(self):
        self.assertEqual(m.midrank([1,1,3,4]),[1.5,1.5,3,4])
        self.assertAlmostEqual(m.spearman([1,2,3],[2,4,6]),1)
        self.assertAlmostEqual(m.spearman([1,2,3],[3,2,1]),-1)
        self.assertIsNone(m.spearman([1,1,1],[2,3,4]))
        with self.assertRaisesRegex(TreeGateError,"NONFINITE_SOURCE_RANK"):
            m.spearman([1,2,float("nan")],[1,2,3])

    def test_within_block_demean_has_declared_singletons_excluded(self):
        r=[{"fold_id":1,"patristic":2,"loss":1},
           {"fold_id":1,"patristic":3,"loss":2},
           {"fold_id":2,"patristic":20,"loss":10}]
        a=m.within_block_corr(r,"patristic","loss")
        self.assertEqual(a["multi_species_blocks_used"],1)
        self.assertEqual(a["species_in_non_singleton_blocks_used"],2)
        self.assertAlmostEqual(a["demeaned_pearson"],1)
        self.assertTrue(a["not_an_independent_cluster_adjusted_regression"])

    def test_duplicate_malformed_and_unverifiable_old_source_stop(self):
        _,codes,data,folds,records=synthetic_records()
        with self.assertRaisesRegex(TreeGateError,"V30_UPSTREAM_CV_SAME_RMSE_METRIC_DRIFT"):
            fake={"endpoints":{
                ep:{"models":{
                    model:{"RMSE_log1p_current_chemical_total":999}
                    for model in MODELS}}
                for ep in m.ENDPOINTS}}
            m.check_existing_model_metrics(records,fake)

    def test_diagnostic_uses_only_frozen_cap_three_and_five(self):
        _,codes,data,folds,records=synthetic_records()
        x=m.summarize_endpoint(records,m.ENDPOINTS[0])
        self.assertEqual(len(x["distance_only_tertiles"]),3)
        self.assertIn("phylo_gap_vs_branch_sq_error",x["rank_correlations_no_significance"])
        self.assertEqual(x["n_original_source_species_scored_once"],31)
        y={"max_clade_3":{"endpoints":{ep:x for ep in m.ENDPOINTS}},
           "max_clade_5":{"endpoints":{ep:x for ep in m.ENDPOINTS}}}
        result=m.diagnostic(y)
        self.assertEqual(len(result["checks"]),4)
        self.assertTrue(result["not_a_causal_distance_effect_or_p_value"])


if __name__=="__main__":
    unittest.main()
