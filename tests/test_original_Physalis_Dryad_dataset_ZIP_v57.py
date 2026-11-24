"""Official *first party* Dryad dataset ZIP is not a guessed source mirror."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_original_Physalis_official_Dryad_dataset_ZIP_v57 as p

class FullOriginalDatasetPackage(unittest.TestCase):
    def test_one_physalis_original_paper_and_version(self):
        self.assertEqual(p.DOI,"10.5061/dryad.xpnvx0kvp")
        _, pin=p.original.check_gate()
        item=pin["source_study_versions"][0]
        self.assertEqual(item["dryad_dataset_id"],176194)
        self.assertEqual(item["published_version"],5)
        self.assertEqual(len(item["author_files"]),10)

    def test_original_2020_source_key_pair_not_2021(self):
        _,pin=p.original.check_gate()
        names={r["file_name"] for r in pin["source_study_versions"][0]["author_files"]}
        self.assertIn("Oikos_2020_cross_season_abundance.csv",names)
        self.assertIn("Oikos_2020_total_plant_yield.csv",names)
        self.assertIn("Oikos_2021_fruit_set_poll_obs.csv",names)

    def test_original_not_proxy_leaking_visits_into_fruit(self):
        self.assertIn("ORIGINAL_PACKAGE_AUTHOR_FILE_DIGEST_UNVERIFIED",
                      open(p.__file__).read())
        self.assertIn("original_reproduction_or_visitation_values_modelled",
                      open(p.__file__).read())

if __name__=="__main__":unittest.main()
