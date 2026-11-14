import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import analyze_target_size_controls_v33 as m
class PairPrecision(unittest.TestCase):
    def test_identical_nonbinary_control_cannot_create_positive_clade(self):
        real={'phylo':np.array([.12345678,.4]),'phylo_genes':np.array([.00765,.1])}
        control={k:np.tile(a,(31,1)) for k,a in real.items()}
        result=m.paired_effect(real,control,[[0],[1]])
        self.assertEqual(result['n_positive_original_clades'],0)
        self.assertEqual(result['n_species_positive_difference'],0)
        self.assertEqual(result['all_original_clade_differences'],[0,0])
if __name__=='__main__':unittest.main()
