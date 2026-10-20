"""Reproduce conditional v5 periodic-capacity results, no field-data input."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from anto_slk_periodic_capacity_v5 import panel


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    parser.add_argument('--expected',type=Path)
    args=parser.parse_args()
    result=panel()
    if args.expected:
        old=json.loads(args.expected.read_text(encoding='utf-8'))
        assert old['schema']==result['schema']
        lookup={(r['phenotype_switching'],r['environment_run_length']):r for r in result['scenarios']}
        for a in old['scenarios']:
            b=lookup[(a['phenotype_switching'],a['environment_run_length'])]
            assert a['phenotype_switching']==b['phenotype_switching']
            assert a['environment_run_length']==b['environment_run_length']
            assert a['regime']==b['regime'] and b['converged']
            assert abs(a['growth_per_generation']-b['growth_per_generation'])<1e-11
            assert abs(a['common_loss_boundary']-b['common_loss_boundary'])<1e-11
            assert max(abs(x-y) for x,y in zip(a['survival_lower'],b['survival_lower']))<1e-10
        for key in ('natural_regain_demonstrated','organism_extinction_inferred',
                    'fixation_or_optimal_strategy_inferred','slk_payoff_calibrated'):
            assert result[key] is False
    text=json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    else: print(text)


if __name__=='__main__': main()
