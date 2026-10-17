"""Reproduce published Cirsium organ contrasts and a separate SYNTHETIC hypothesis.

Never treat the 2020 organ-expression study as measured Ryukyu competence.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_cirsium_organ_evidence_v6 import report
from anto_modular_competence_v6 import synthetic_organ_protection_grid


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,
                        default=ROOT/"data/molecular/cirsium_japonicum_organ_evidence_v6.json")
    parser.add_argument("--output",type=Path,default=None)
    args=parser.parse_args()
    result={
        "research_version":"anto_v6",
        "measured_organ_evidence_from_published_study":report(args.source),
        "theory_hypothesis_not_fitted_to_source":synthetic_organ_protection_grid(),
        "hard_separation":{
            "organ_evidence_is_empirical":True,
            "illustrative_transition_rates_are_empirical":False,
            "Ryukyu_floral_reexpression_confirmed":False,
            "published_2020_Cirsium_is_independent_Ryukyu_reversion":False,
            "source_study_natural_clade_replicates":1,
            "no_universal_law_claim":True
        }
    }
    text=json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding="utf-8")
    else: print(text)
    return result


if __name__=="__main__":
    main()
