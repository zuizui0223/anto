"""Exact model-equivalence and repeat-expression reporting with no plant data."""
from pathlib import Path
import json
import argparse
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_afterlife_induction_identifiability_v10b import panel

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    data=panel()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("INDUCTION_IDENTIFIABILITY_PASS",
        json.dumps({"same_initial_positive":data["first_induction_counterexample"]["first_pulse_visible_positive_both"],
                    "independent_genetic_module_assay_needed":True,
                    "longitudinal_persistence_needed":True,
                    "natural_trait_regain_demonstrated":False}))

if __name__=="__main__":main()
