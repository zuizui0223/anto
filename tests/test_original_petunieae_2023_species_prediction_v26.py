"""Data identity, denominator, no leakage and deterministic CV tests, v26."""
import io
import math
import sys
import unittest
from copy import deepcopy
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_original_petunieae_2023_species_prediction_v26 import (
    ANTHOCYANIDINS,FLAVONOLS,MODELS,RIDGE_LAMBDA,OUTGROUP,
    parse_named_csv,numeric,reconcile_original,
    linear_solve,ridge_predict,blocked_folds,evaluate
)


def synthetic_original_rows():
    """Fictional 59 ingroups + single BROW, exact 3-source row key match."""
    types=[f"AA{i:02d}" for i in range(20)]+[f"BB{i:02d}" for i in range(19)]+[
        f"CC{i:02d}" for i in range(20)]+["BROW"]
    chemical=list(ANTHOCYANIDINS+FLAVONOLS)
    genes=sorted(set().union(*map(set,MODELS.values())))
    all_header=["key_0"]+genes+chemical
    gene_header=[""]+genes+["ODO1"]
    pigment_header=["species"]+chemical
    rows=[];gene_rows=[];chemical_rows=[]
    for i,code in enumerate(types):
        values=[str(1+i) for _ in genes]
        pigments=[str((i+1)*j/100) for j in range(1,len(chemical)+1)]
        rows.append([code]+values+pigments)
        gene_rows.append([code]+values+["1"])
        chemical_rows.append([code]+pigments)
    def csv(h,data):
        return ("\n".join([",".join(h)]+[",".join(r) for r in data])+"\n").encode()
    return csv(all_header,rows),csv(gene_header,gene_rows),csv(pigment_header,chemical_rows)


def sample_panel():
    """15 fake taxa with three distinct 2-letter folds and exact relation."""
    result={}
    for g in ("AA","BB","CC"):
        for i in range(5):
            k=g+str(i)
            x=(i+1)+({"AA":0,"BB":3,"CC":6}[g])
            v=math.log1p(x)
            result[k]={"genes":{"CHI":v,"CHS":v/2,"F3H":v/3,
                       "DFR":v/4,"ANS":v/5,"F3primeH":v/6,
                       "F35H":v/7,"FLS":v/8},
                       "responses":{"anthocyanidins_log1p_total":2*v,
                                    "flavonols_log1p_total":3*v}}
    return result


class OriginalPetunieaeProcessedSpeciesHoldoutV26(unittest.TestCase):
    def test_original_2023_input_endpoint_contract(self):
        self.assertEqual(OUTGROUP,"BROW")
        self.assertEqual(len(ANTHOCYANIDINS),6)
        self.assertEqual(len(FLAVONOLS),3)
        self.assertEqual(len(MODELS),6)
        self.assertEqual(RIDGE_LAMBDA,10.0)
        self.assertNotIn("BROW",list(ANTHOCYANIDINS)+list(FLAVONOLS))
        self.assertTrue(all(not any("_mgg" in z for z in genes)
                            for genes in MODELS.values()))

    def test_fixed_59_species_matched_processed_files_and_outgroup_excluded(self):
        raw,genes,pig=synthetic_original_rows()
        source,info=reconcile_original(raw,genes,pig)
        self.assertEqual(len(source),59)
        self.assertNotIn("BROW",source)
        self.assertTrue(info["Browallia_americana_outgroup_removed"])
        self.assertEqual(info["source_original_pigment_cells_verified_by_independent_table"],540)
        self.assertEqual(info["source_original_gene_cells_verified_by_independent_table"],480)
        self.assertTrue(info["no_visible_hue_observed_baseline_in_this_table"])

    def test_mislabeled_original_source_species_row_is_rejected(self):
        raw,genes,pig=synthetic_original_rows()
        with self.assertRaisesRegex(ValueError,"do not match"):
            reconcile_original(raw,genes.replace(b"AA00",b"XXXX"),pig)
        with self.assertRaisesRegex(ValueError,"do not match"):
            reconcile_original(raw,genes,pig.replace(b"BB00",b"XXXX"))

    def test_original_pigment_source_reconciliation_fails_on_mismatched_cell(self):
        raw,genes,pig=synthetic_original_rows()
        # Change exactly the first chemical value, while preserving shape.
        lines=pig.decode().splitlines()
        row=lines[1].split(",");row[1]="9.999";lines[1]=",".join(row)
        with self.assertRaisesRegex(ValueError,"HPLC rows disagree"):
            reconcile_original(raw,genes,("\n".join(lines)+"\n").encode())

    def test_original_gene_source_reconciliation_fails_on_mismatched_cell(self):
        raw,genes,pig=synthetic_original_rows()
        lines=genes.decode().splitlines()
        row=lines[1].split(",");row[1]="9999.99";lines[1]=",".join(row)
        with self.assertRaisesRegex(ValueError,"gene rows disagree"):
            reconcile_original(raw,("\n".join(lines)+"\n").encode(),pig)

    def test_original_source_missing_zero_and_nan_are_not_silently_imputed(self):
        row={"CHI":"","DFR":"NaN","ANS":"-1","CHS":"0"}
        for v in ("CHI","DFR","ANS"):
            with self.assertRaises(ValueError):
                numeric(row,v)
        self.assertEqual(numeric(row,"CHS"),0.)

    def test_source_rectangular_60_species_and_column_identity(self):
        raw,_,_=synthetic_original_rows()
        one,h,key=parse_named_csv(raw,60)
        self.assertEqual(len(one),60)
        self.assertEqual(key,"key_0")
        self.assertIn("Del_mgg",h)
        with self.assertRaisesRegex(ValueError,"row count"):
            parse_named_csv(raw,59)
        with self.assertRaisesRegex(ValueError,"column names"):
            parse_named_csv(b"species,CHI,CHI\nA,1,2\n",1)

    def test_independent_small_dense_linear_solver_and_noninvertible_matrix(self):
        self.assertEqual(linear_solve([[2,0],[0,3]],[2,6]),[1.,2.])
        with self.assertRaisesRegex(ValueError,"Singular"):
            linear_solve([[0,0],[0,0]],[1,1])
        self.assertEqual(linear_solve([],[]),[])

    def test_training_only_intercept_not_whole_original_test_response(self):
        s=sample_panel()
        test=s["AA0"]
        train={k:v for k,v in s.items() if k!="AA0"}
        a=ridge_predict(train,test,"anthocyanidins_log1p_total",())
        expected=sum(v["responses"]["anthocyanidins_log1p_total"]
                     for v in train.values())/len(train)
        self.assertAlmostEqual(a,expected)
        test_mutant=deepcopy(test)
        test_mutant["responses"]["anthocyanidins_log1p_total"]=-1e7
        self.assertEqual(a,ridge_predict(
            train,test_mutant,"anthocyanidins_log1p_total",()))
        self.assertEqual(ridge_predict(train,test,"flavonols_log1p_total",()),
                         sum(z["responses"]["flavonols_log1p_total"] for z in train.values())/len(train))

    def test_training_only_standardized_ridge_and_heldout_target_not_used(self):
        samples=sample_panel()
        test=samples["CC4"]
        train={k:v for k,v in samples.items() if k!="CC4"}
        predicted=ridge_predict(train,test,"anthocyanidins_log1p_total",("CHI",))
        self.assertTrue(math.isfinite(predicted))
        modified=deepcopy(test)
        modified["responses"]["anthocyanidins_log1p_total"]=1e8
        self.assertEqual(predicted,ridge_predict(
            train,modified,"anthocyanidins_log1p_total",("CHI",)))
        self.assertNotEqual(predicted,ridge_predict(
            train,test,"anthocyanidins_log1p_total",()))
        self.assertNotEqual(predicted,test["responses"]["anthocyanidins_log1p_total"])

    def test_species_and_prefix_group_folds_never_reuse_test_species_as_train(self):
        samples=sample_panel()
        loo=blocked_folds(samples,"species_LOO")
        grouped=blocked_folds(samples,"first_two_source_code_characters_NOT_accepted_genus")
        self.assertEqual(len(loo),15)
        self.assertEqual([len(z) for z in grouped],[5,5,5])
        for folds in (loo,grouped):
            self.assertEqual({x for fold in folds for x in fold},set(samples))
            self.assertEqual(len([x for fold in folds for x in fold]),len(samples))
        with self.assertRaisesRegex(ValueError,"Unregistered"):
            blocked_folds(samples,"full_phylogenetic_clade_LOO")

    def test_one_original_source_clade_no_species_recovery_CV_claim(self):
        samples=sample_panel()
        res=evaluate(samples,"species_LOO")
        self.assertEqual(res["n_folds"],15)
        self.assertTrue(res["refitted_train_only_scaling_and_intercept"])
        self.assertTrue(res["code_prefix_is_NOT_verified_taxonomic_genus_or_an_independent_regain_clade"])
        for label,row in res["endpoints"].items():
            self.assertEqual(row["fixed_unretuned_ridge_lambda"],10.0)
            self.assertEqual(row["models"]["train_mean_baseline"]["relative_RMSE_to_source_LOO_baseline"],1.)
            self.assertLess(row["models"]["CHI_only"]["out_of_fold_RMSE_log1p_mg_per_g"],
                            row["models"]["train_mean_baseline"]["out_of_fold_RMSE_log1p_mg_per_g"])

    def test_prefix_group_holdout_fails_if_source_grouping_degenerate(self):
        samples={f"AA{i}":data for i,data in enumerate(sample_panel().values())}
        with self.assertRaisesRegex(ValueError,"grouping"):
            blocked_folds(samples,"first_two_source_code_characters_NOT_accepted_genus")


if __name__=="__main__":
    unittest.main()
