import sys,json,itertools,random,argparse
import numpy as np
from fractions import Fraction as F
from pathlib import Path
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path)
args=parser.parse_args()
sys.path.insert(0,str(root/'src'))
from anto_slk_recurrent_capacity_v4 import CapacityProcess
from anto_slk_periodic_capacity_v5 import PeriodicCapacity
errs=[]; agreements=0;near=0
periods=[]
for high,low,switch,loss,k in itertools.product(['.75','.85','.95'],['.25','.4','.6'],['.1','.5','.9'],[0,'.1','.2','.3'],[1,2,4]):
 a=CapacityProcess(high,low,loss,switch,switch)
 b=CapacityProcess(low,high,loss,switch,switch)
 periods.append(PeriodicCapacity((a,)*k+(b,)*k))
rng=random.Random(20261020)
for i in range(40):
 phases=tuple(CapacityProcess(F(rng.randint(3,9),10),F(rng.randint(3,9),10),F(rng.randint(0,3),10),F(rng.randint(1,9),10),F(rng.randint(1,9),10)) for _ in range(rng.randint(2,6)))
 periods.append(PeriodicCapacity(phases))
for p in periods:
 q=np.eye(2)
 for phase in p.phases: q=q@np.array(phase.mean_matrix(),dtype=float)
 cycle=max(np.linalg.eigvals(q).real)
 growth=cycle**(1/len(p.phases))
 errs.append(abs(growth-p.growth_per_generation()))
 if abs(cycle-1)<1e-10: near+=1
 else:
  expected='supercritical' if cycle>1 else 'subcritical'
  assert p.regime()==expected,(p,cycle,p.regime())
  agreements+=1
assert max(errs)<1e-12
sample=[p for p in periods if p.regime()=='supercritical' and p.growth_per_generation()>1.005][::15]
checks=[]
for p in sample:
 s=p.survival()
 assert s['converged'] and s['width']<1e-11 and s['residual']<1e-11
 checks.append(s['iterations'])
result={'grid_cases':len(periods),'symmetric_grid_cases':324,'asymmetric_random_cases':40,'numpy_eigenvalue_classification_agreements':agreements,'near_critical_float_comparisons_skipped':near,'max_per_generation_absolute_difference':max(errs),'sampled_supercritical_survival_checks':len(checks),'max_iterations_in_sample':max(checks),'input_base_commit':'0d57f461c3e2e00181251c412f32d2a4df9c898e'}
if args.output:
 args.output.parent.mkdir(parents=True,exist_ok=True)
 args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
