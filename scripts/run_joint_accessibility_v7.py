"""Reproduce the exact resource-inventory Fréchet counterexamples and migration timing.

Only hypothetical state probabilities/rates. No real Cirsium phenotype,
flower-stage RNA, genomic ancestral origin, or cross-angiosperm biological
generality is inferred here. All empirical Cirsium design lives in aza3.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_joint_accessibility_v7 import exact_report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=None)
    args=parser.parse_args()
    data=exact_report()
    data["program_role"]="generic_trait_biogeography_theory"
    data["Cirsium_empirical_owner"]="zuizui0223/aza3"
    data["validated_in_vivo_colonization_effect"]=False
    data["validated_capitulum_LEGO_recombination"]=False
    output=json.dumps(data,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output,encoding="utf-8")
    else:
        print(output)
    return data


if __name__=="__main__":
    main()
