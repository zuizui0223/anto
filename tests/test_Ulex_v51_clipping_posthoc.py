"""Post-v51-results original Ulex clipping saturation controls."""
import json
import os
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_original_Ulex_v51_clipping_posthoc as c
import analyze_original_Ulex_elevation_analogues_v51 as v


class ClipPosthoc(unittest.TestCase):
    def test_disclosed_post_outcome_and_parent_blob(self):
        r=c.lock()
        self.assertTrue(r["retrospective_posthoc"])
        self.assertFalse(r["valid_randomization_p_value"])
        self.assertEqual(c.git_blob_sha1(c.PARENT.read_bytes()),
                         "8f1d5435d2d14c2b2c999671ff9801cffab01a71")

    def test_six_original_sites_not_six_genetically_independent_pops(self):
        p=c.read_parent(c.lock())
        self.assertTrue(p["six_sites_within_one_continuous_original_population"])
        self.assertEqual(p["n_original_site_holdouts"],6)
        self.assertEqual(p["n_exact_train_two_site_pair_fits"],60)

    def test_preclip_prediction_without_artificial_truncation(self):
        fit={"beta":[.5,1.0],"stats":{
            "elev":(0.,1.),"area":(0.,1.),"weight":(0.,1.)},
            "traits":False}
        row={"elev":3.0,"area":1.0,"weight":1.0}
        self.assertAlmostEqual(c.preclip(fit,row),3.5)
        row["elev"]=-5.0
        self.assertAlmostEqual(c.preclip(fit,row),-4.5)

    def test_source_record_disallows_biomolecular_or_visitor_inference(self):
        r=c.lock()
        self.assertEqual(r["independently_sampled_original_studies_added"],0)
        self.assertFalse(r["external_pollinator_individual_joined"])
        self.assertEqual(r["models"],["M0","M1"])

    def test_six_site_recomputation_if_original_md5_staged(self):
        root=os.environ.get("ULEX_V51_SOURCE_ROOT")
        if not root:
            self.skipTest("Original Zenodo bytes staged only under CI")
        reg=v.source_registry()
        raw=v.download_source(reg,source_root=root)
        source,audit=v.parse_original(raw,reg)
        out=c.source_audit(source,c.read_parent(c.lock()),c.lock())
        self.assertEqual(len(out["derived_original_site_diagnostics"]),18)
        self.assertEqual(set(out["aggregated_all_sites"]),set(c.REGIMES))
        for entry in out["aggregated_all_sites"].values():
            self.assertEqual(entry["n_original_heldout_individuals_once_each"],202)
            for m in c.MODELS:
                q=entry["by_model"][m]
                self.assertEqual(q["n_raw_predictions_below_zero"]+
                                 q["n_raw_predictions_above_one"]+
                                 q["n_raw_predictions_within_unit"],202)

    def test_source_full_posthoc_receipt_recomputed_after_freeze(self):
        root=os.environ.get("ULEX_V51_SOURCE_ROOT")
        if not root or not c.RECEIPT.exists():
            self.skipTest("Exact source bytes and frozen source receipt required")
        reg=v.source_registry()
        data,audit=v.parse_original(v.download_source(reg,source_root=root),reg)
        actual=c.source_audit(data,c.read_parent(c.lock()),c.lock())
        self.assertEqual(actual,json.loads(c.RECEIPT.read_text()))
        self.assertTrue(actual["posthoc_not_primary_confirmation"])


if __name__=="__main__":
    unittest.main()
