"""Reproduce preregistered-style synthetic checks; NEVER interpret as field data."""
import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from anto_colonization import Parameters, Plant, simulate, first_pigmented, first_pigmented_route, order_experiment_probability
from anto_identifiability import benchmark


def wilson(successes, total, z=1.959963984540054):
    if total <= 0:
        raise ValueError('Total must be positive')
    p = successes / total
    a = z*z/total
    c = (p+a/2)/(1+a)
    h = z*math.sqrt((p*(1-p)+a/4)/total)/(1+a)
    return [max(0., c-h), min(1., c+h)]


def scenarios():
    return {
        'competent_white_island': (
            Parameters(('island',), {'island':30}, fecundity=2.5,
                       reexpression=.08, degradation={'island':.06}),
            {'island':[Plant('local_white','W+')]*12}),
        'degraded_white_island': (
            Parameters(('island',), {'island':30}, fecundity=2.5,
                       reexpression=.08, degradation={'island':.06}),
            {'island':[Plant('local_white','W-')]*12}),
        'pigmented_seed_immigration': (
            Parameters(('source','island'), {'source':30,'island':30},
                       seed_migration={('source','island'):.06},
                       fecundity=2.5, reexpression=0),
            {'source':[Plant('immigrant','P')]*12,
             'island':[Plant('local_white','W-')]*12}),
        'pollen_introgression_only': (
            Parameters(('source','island'), {'source':30,'island':30},
                       pollen_delivery={('source','island'):.5},
                       fecundity=2.5, reexpression=0,
                       introgression=.15),
            {'source':[Plant('donor','P')]*12,
             'island':[Plant('local_white','W-')]*12}),
    }


def main(replicates=400, horizon=8):
    if replicates <= 0 or horizon <= 0:
        raise ValueError('Replicates and horizon must be positive')
    output = {
        'version':'synthetic_benchmark_v2', 'seed_base': 20261017,
        'empirical_data_admitted': False,
        'interpretation':'Toy-model consequences under stipulated rates, NOT Cirsium rates',
        'replicates_per_scenario':replicates, 'generations':horizon,
        'scenarios':{},
        'order_test':{
            'parameters':{'delta_source':.02,'delta_island':.2,
                          'horizon':8,'reexpression':.1,'founder_count':8},
            'by_delay':{str(d):order_experiment_probability(.02,.2,d,8,.1,8)
                        for d in [0,5,10,20]}},
        'history_identifiability':benchmark(seed=20261017),
    }
    for key, (par, initial) in scenarios().items():
        total_coloured_end = 0
        total_coloured_ever = 0
        event_totals = {}
        for i in range(replicates):
            result = simulate(par, initial, horizon, seed=20261017+i)
            total_coloured_end += int(result['histories'][-1]['regions']['island']['pigmented']>0)
            total_coloured_ever += int(first_pigmented(result,'island') is not None)
            for event in result['events']:
                for route, count in event['regions']['island'].items():
                    event_totals[route] = event_totals.get(route,0)+count
        output['scenarios'][key] = {
            'final_coloured_fraction':total_coloured_end/replicates,
            'final_coloured_wilson95':wilson(total_coloured_end,replicates),
            'ever_coloured_fraction':total_coloured_ever/replicates,
            'ever_coloured_wilson95':wilson(total_coloured_ever,replicates),
            'accepted_recruitment_event_counts':dict(sorted(event_totals.items()))}
    mixed = Parameters(('source','island'), {'source':30,'island':30},
                       seed_migration={('source','island'):.005},
                       pollen_delivery={('source','island'):.4},
                       fecundity=1.6, reexpression=.005,
                       degradation={'island':.08}, introgression=.005)
    source = {'source':[Plant('source_P','P')]*8,
              'island':[Plant('island_W','W+')]*8}
    route_counts = {}
    for i in range(replicates):
        simulation = simulate(mixed, source, horizon, seed=20261017+i)
        route = first_pigmented_route(simulation, 'island')['route']
        route_counts[route] = route_counts.get(route, 0) + 1
    output['mixed_routes'] = {
        'note':'First coloured island cohort: outcomes from a mixture of all three mechanisms; same-generation ties retained',
        'first_route_counts':dict(sorted(route_counts.items())),
        'first_route_proportions':{k:v/replicates for k,v in sorted(route_counts.items())},
    }
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--replicates',type=int,default=400)
    parser.add_argument('--horizon',type=int,default=8)
    parser.add_argument('--output',type=Path,default=None)
    args = parser.parse_args()
    results = main(args.replicates,args.horizon)
    content = json.dumps(results,ensure_ascii=False,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content, encoding='utf-8')
    else:
        print(content)
