"""Regression tests: original Mertensia source identity and entire-domain folds."""
import math
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_Mertensia_season_water_transfer_v59 as v

class ContextTransfer(unittest.TestCase):
    def test_pinned_v58_source(self):
        r=v.registry()
        self.assertEqual(r["n_plants"],105)
        self.assertEqual(r["v58_result_git_blob"],"cd2de70901ba305c05c86ccc0aa4dd08d803014d")
        self.assertEqual(r["new_original_studies"],0)
        self.assertFalse(r["causal_claim"])
    @staticmethod
    def fixture():
        rows=[]
        for w in ("1","2","3","4"):
            for a in ("1","2","3","4","5"):
                for j in range(6 if w=="1" else 5):
                    i=len(rows)
                    t=v.TREAT[(i+j)%3]
                    rows.append({"id":(w,a,t,str(i)),"week":w,"array":a,
                                 "group":(w,a),"treat":t,"week_num":float(w),
                                 "log_flower":math.log1p(i%17+1),
                                 "mcorW":2.+i*.02,"mcorL":5.+(i%13)*.13,
                                 "log_visit":math.log1p(i%9),
                                 "mVWC":10.+i%18,"y":math.log1p(i%4)})
        assert len(rows)==105
        return rows
    def test_all_policies_exact_partition(self):
        rows=self.fixture()
        for policy,count in (("array20",20),("safe16",16),("week4",4),("water3",3)):
            folds=v.splits(rows,policy)
            self.assertEqual(len(folds),count)
            self.assertEqual(sorted(i for _,_,te in folds for i in te),list(range(105)))
    def test_week4_heldout_together(self):
        rows=self.fixture()
        k,tr,te=next(x for x in v.splits(rows,"safe16") if x[0]=="4|ALL")
        self.assertTrue(all(rows[i]["week"]!="4" for i in tr))
        self.assertTrue(all(rows[i]["week"]=="4" for i in te))
    def test_no_unseen_water_dummy_or_week_dummy(self):
        rows=self.fixture()
        f=v.fit(rows,"linear_week_vwc",True,"water3")
        self.assertEqual(len(f["beta"]),1+len(f["fs"]))
        f=v.fit(rows,"linear_week",True,"week4")
        self.assertEqual(len(f["beta"]),3+len(f["fs"]))
        self.assertIn("week_num",f["fs"])
        self.assertNotIn("mVWC",f["fs"])
        self.assertNotIn("y",f["fs"])
    def test_gain_decomposition(self):
        x=v.errors({"y":3},1,2)
        z=v.agg([x])
        self.assertEqual(z["gain"],3)
        self.assertEqual(z["alignment_2ed"],4)
        self.assertEqual(z["correction_d2"],1)
        self.assertEqual(v.agg([v.errors({"y":3},1,6)])["gain"],-5)
    def test_target_outcome_not_used_in_prediction(self):
        rows=self.fixture()
        f=v.fit(rows[1:],"linear_week_vwc",True,"week4")
        x=dict(rows[0])
        p=v.predict(f,x)
        x["y"]=999.
        self.assertEqual(p,v.predict(f,x))
    def test_cross_week_physical_plant_id_stop(self):
        rows=self.fixture()
        rows[-1]["id"]=tuple(rows[-1]["id"][:3])+(rows[0]["id"][3],)
        with self.assertRaisesRegex(v.V59Error,"POTENTIAL_SAME_PHYSICAL_PLANT_LEAKAGE"):
            v.splits(rows,"week4")
    def test_equal_week_and_plant_weight_reported(self):
        rows=self.fixture()
        a=v.evaluate(rows,"week4","linear_week")
        self.assertEqual(a["n_plants"],105)
        self.assertEqual(a["n_folds"],4)
        self.assertEqual(set(a["by_week"]),set(v.WEEKS))
        self.assertAlmostEqual(a["fold_equal"]["gain"],a["week_equal"]["gain"])
        self.assertAlmostEqual(a["plant_weighted"]["gain"],
            a["plant_weighted"]["alignment_2ed"]-a["plant_weighted"]["correction_d2"])
    def test_negative_training_sd_fails(self):
        rows=self.fixture()
        for x in rows:x["mVWC"]=20.
        with self.assertRaisesRegex(v.V59Error,"ZERO_TRAIN_VARIANCE_mVWC"):
            v.fit(rows,"linear_week_vwc",True,"week4")
    def test_200_placebos_preserve_registered_group_grain(self):
        rows=self.fixture()
        actual=v.evaluate(rows,"week4","linear_week")
        p=v.placebo(rows,actual,v.registry())
        self.assertEqual(p["n"],200)
        self.assertTrue(p["reference_is_not_a_valid_causal_p_value"])
if __name__=="__main__":unittest.main()
