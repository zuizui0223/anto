"""Recompute route-specific recovery windows; no empirical rates are fitted."""
import argparse
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import anto_afterlife_route_window_v12 as model


def report():
    loss=model.Loss('.02','.30','.01')
    rates=model.Rates('.40','.20',loss)
    peak=model.alternative_peak(rates,loss)
    ages=(0.,5.,20.,peak['age'])
    cases=[]
    for age in ages:
        cases.append({'silent_age_model_units':age,
                      'eventual':model.after_delay(rates,loss,age),
                      'followup_10_model_units':model.after_delay(rates,loss,age,10.),
                      'alternative_ever_no_first_event_censoring':
                          model.alternative_ever_without_preemption(rates,loss,age)})
    grid=model.parameter_grid()
    interior=[r for r in grid if r['peak_kind']=='interior']
    symmetry=model.Loss(.02,.1,.1)
    near=[]
    for factor in (.9,1.,1.1):
        r=model.Rates(.4*factor,.2,loss)
        near.append({'factor_on_ancestral_activation':factor,
                     'peak':model.alternative_peak(r,loss)})
    source=ROOT/'src/anto_afterlife_route_window_v12.py'
    output={
        'schema':'anto_afterlife_route_window_v12',
        'analysis_date':'2026-10-20',
        'base_anto_commit':'3efa90dbd5f95469ff4f8732606e8c8d64b0a347',
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'synthetic_only':True,
        'parameters_estimated_from_real_plants':False,
        'endpoint':'Identity and probability of FIRST re-expression after opportunity returns',
        'absorbing_labels_mean_first_event_record_not_permanent_biological_state':True,
        'initial_condition':'All C,A,B functional resources intact; expressed phenotype absent',
        'pre_opportunity':'No re-expression, independent irreversible C/A/B losses',
        'parameter_units':'Arbitrary model time; not calibrated years or generations',
        'rates':{'activation_ancestral':.4,'activation_alternative':.2,
                 'loss_common':.02,'loss_ancestral':.30,'loss_alternative':.01,
                 'same_loss_rates_before_and_after_opportunity':True},
        'no_new_module_creation_mutation_import_or_regain_of_lost_resource':True,
        'peak_for_eventual_alternative_first':peak,
        'reference_cases':cases,
        'strict_local_neighborhood':near,
        'deterministic_parameter_grid':{
            'settings':len(grid),'interior_maxima':len(interior),
            'all_interior_have_higher_absolute_alternative_first':all(r['alternative_at_peak']>r['alternative_at_zero'] for r in interior),
            'all_interior_have_nonincreasing_total_first_recovery':all(r['any_at_peak']<=r['any_at_zero']+1e-12 for r in interior),
            'grid_is_not_a_probability_distribution_or_natural_incidence_estimate':True},
        'negative_control':{
            'symmetric_routes_peak':model.alternative_peak(model.Rates(.2,.2,symmetry),symmetry),
            'alternative_ever_counted_after_ancestral_removes_window':True,
            'interpreting_B_first_increase_as_faster_B_activation_is_invalid':True},
        'biological_establishment_fixation_fitness_and_natural_regain_measured':False,
        'novelty_boundary':'Standard competing-risk mathematics; new route-resolved conditional extension in this repository, not a new general probability theorem',
        'empirical_owner_Cirsium':'zuizui0223/aza3',
        'external_motivation_not_calibration':[
            'https://doi.org/10.1093/plcell/koab114',
            'https://doi.org/10.1093/molbev/msy117'],
        'reproducibility':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
    }
    return output,grid


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--grid-csv',type=Path)
    args=p.parse_args()
    output,grid=report()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    if args.grid_csv:
        args.grid_csv.parent.mkdir(parents=True,exist_ok=True)
        with args.grid_csv.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=grid[0].keys())
            writer.writeheader();writer.writerows(grid)
    print(json.dumps({'model_only':True,'peak':output['peak_for_eventual_alternative_first'],
                      'grid':output['deterministic_parameter_grid']},sort_keys=True))

if __name__=='__main__': main()
