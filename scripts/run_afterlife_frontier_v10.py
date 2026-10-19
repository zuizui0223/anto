"""Produce a deterministic, exact-fraction generic afterlife frontier receipt."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_afterlife_frontier_v10 import panel


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path)
    args=a.parse_args()
    result=panel()
    result["empirical_generalization_demonstrated"]=False
    result["Cirsium_empirical_owner"]="zuizui0223/aza3"
    content=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding="utf-8")
    else:
        print(content)
    return result


if __name__=="__main__":main()
