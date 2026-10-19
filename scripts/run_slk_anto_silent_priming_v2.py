"""Write the exact ANTO-SLK silent priming opportunity result as JSON.

No estimated biological rates; all fractions represent synthetic probabilities.
"""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_slk_silent_priming_window_v2 import panel


def normalize(value):
    if isinstance(value,Fraction):
        return str(value)
    if isinstance(value,dict):
        return {str(k):normalize(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [normalize(v) for v in value]
    return value


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    receipt=normalize(panel())
    content=json.dumps(receipt,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding="utf-8")
    else:
        print(content)
    return receipt


if __name__=="__main__":
    main()
