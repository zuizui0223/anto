"""Run source-admission check plus exact *hypothetical* observational designs."""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from anto_assay_design import run_design_grid
from anto_four_history_design import comparison
from anto_evidence_gate import report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--registry',type=Path,default=ROOT/'data'/'empirical_sources_v3.json')
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args()
    result={
       'version':'anto_empirical_gate_and_assay_design_v3',
       'date':'2026-10-17',
       'empirical_inference_performed':False,
       'literature_source_audit':report(args.registry),
       'three_history_designs':run_design_grid(),
       'four_history_designs':comparison(),
       'claim_boundary':'The rank computations are exact given hypothetical known independent emissions; not an empirical fitted history or validated molecular assay.'
    }
    s=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(s,encoding='utf-8')
    else:
        print(s)
    return result

if __name__=='__main__':
    main()
