"""Emit the exact rational, theory-only SLK × ANTO ecological window receipt."""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_slk_ecological_window_v1 import panel


def serializable(value):
    if isinstance(value,Fraction):
        return str(value)
    if isinstance(value,dict):
        return {k:serializable(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [serializable(v) for v in value]
    return value


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=serializable(panel())
    content=json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding="utf-8")
    else:
        print(content)
    return result


if __name__=="__main__":
    main()
