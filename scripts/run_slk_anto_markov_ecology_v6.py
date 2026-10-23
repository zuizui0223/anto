"""Emit deterministic, sampled-and-finite stochastic ecology receipts.

No sequence data, flora or empirically fitted environmental rates are used.
All modeled events concern competent descendants, not trait fixation.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_slk_markov_ecology_v6 import panel


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    p=panel()
    payload=json.dumps(p,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(payload,encoding="utf-8")
    else:
        print(payload)
    return p


if __name__=="__main__":
    main()
