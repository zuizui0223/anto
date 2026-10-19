"""Reproduce anto v9 conditional historical-memory experiment.

Pure synthetic exact probabilities; no empirical floral re-expression rates.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_history_state_sufficiency_v9 import comparative_sensitivity_grid


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    result=comparative_sensitivity_grid()
    result["Cirsium_empirical_owner"]="zuizui0223/aza3"
    result["experimental_success_rate_calibrated_in_plants"]=False
    content=json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(exist_ok=True,parents=True)
        args.output.write_text(content,encoding="utf-8")
    else:
        print(content)
    return result


if __name__=="__main__":
    main()
