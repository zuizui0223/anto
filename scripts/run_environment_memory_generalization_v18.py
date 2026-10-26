"""Reproduce exact non-symmetric shared-environment branching contrasts.

All numbers are constructed fractions, not species-specific estimates.
"""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from anto_environment_memory_generalization_v18 import panel


def serializable(x):
    if isinstance(x, Fraction):
        return str(x)
    if isinstance(x, dict):
        return {str(k): serializable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [serializable(i) for i in x]
    return x


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args=parser.parse_args()
    result=serializable(panel())
    content=json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding='utf-8')
    else:
        print(content)


if __name__=='__main__':
    main()
