"""Frozen exact-rational v17 receipt; no empirical parameter fitting."""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_environment_horizon_bridge_v17 import bridge_panel


def serialize(x):
    if isinstance(x,Fraction):
        return str(x)
    if isinstance(x,dict):
        return {str(k):serialize(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [serialize(y) for y in x]
    return x


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    data=json.dumps(serialize(bridge_panel()),ensure_ascii=False,
                    sort_keys=True,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(data,encoding="utf-8")
    else:
        print(data)


if __name__=="__main__":
    main()
