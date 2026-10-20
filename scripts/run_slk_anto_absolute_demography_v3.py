"""Export a deterministic ANTO × SLK absolute-demography comparison.

The PGF coefficient recurrence is finite and exact in form; numerical
values are evaluated with 65-digit Decimal precision, not large fractions.
This is a constructed model, not a measurement of recovery or rescue.
"""
import argparse
import json
import sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_slk_absolute_demography_v3 import panel


def encode(obj):
    if isinstance(obj,Fraction):
        return str(obj)
    if isinstance(obj,Decimal):
        return format(obj, ".18f")
    if isinstance(obj,dict):
        return {str(k):encode(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)):
        return [encode(v) for v in obj]
    return obj


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    content=json.dumps(encode(panel()),sort_keys=True,ensure_ascii=False,
                       indent=2)+"\n"
    if args.output is None:
        print(content)
    else:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding="utf-8")


if __name__=="__main__":
    main()
