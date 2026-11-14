"""ANTO v31: authenticated species-LOO pigment prediction and distance audit.

One inspected 2023 chemical study; no causal distance effect or new biological
replication. Parent modules are imported only in report(), so numerical and
leakage controls run without network access or real chemical observations.
"""
from __future__ import annotations
import argparse
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
from statistics import median

TOLERANCES = (1e-6, 1e-5, 1e-4)
PRIMARY_TOL = 1e-5
BRANCH_GENES = ('F3primeH', 'F35H', 'FLS')
ENDPOINTS = ('anthocyanidins_log1p_total', 'flavonols_log1p_total')
PARENT = '967dfb39c048301612f291d102309720f9ba31e0'
RECEIPTS = {
    '26': 'petunieae_2023_original_species_pigment_prediction_v26_20261028.json',
    '27': 'petunieae_original_tree_clade_holdout_v27_20261114.json',
    '28': 'petunieae_alternative_published_phylogenies_v28_20261114.json',
    '29': 'petunieae_fine_discordance_label_reference_v29_20261114.json',
    '30': 'petunieae_patristic_prediction_gaps_v30_20261114.json',
}


def tk(tol):
    return f'{tol:.0e}'


def finite(value, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError('NONNUMERIC_VALUE')
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError('NONFINITE_OR_NEGATIVE_VALUE')
    return value


def tree_distances(root):
    """Independent Decimal path sums; no pruning and no root-stem edge."""
    paths = {}
    seen = set()
    def visit(node, ancestors):
        if id(node) in seen:
            raise ValueError('TREE_REUSED_NODE_OR_CYCLE')
        seen.add(id(node))
        if node.children:
            for child in node.children:
                finite(child.length, nonnegative=True)
                visit(child, ancestors + [(id(child), Decimal(str(child.length)))])
        else:
            if not node.label or node.label in paths:
                raise ValueError('DUPLICATE_OR_EMPTY_TIP')
            paths[node.label] = dict(ancestors)
    visit(root, [])
    if len(paths) < 2:
        raise ValueError('TOO_FEW_TIPS')
    result = {}
    for a, path_a in paths.items():
        result[a] = {}
        for b, path_b in paths.items():
            unique_a = set(path_a) - set(path_b)
            unique_b = set(path_b) - set(path_a)
            distance = sum((path_a[k] for k in unique_a), Decimal(0))
            distance += sum((path_b[k] for k in unique_b), Decimal(0))
            result[a][b] = float(distance)
    if not any(d > 0 for row in result.values() for d in row.values()):
        raise ValueError('ALL_SOURCE_DISTANCES_ZERO')
    return result


def distance_bands(gaps, tolerance):
    """Complete-link sorted bands, never chaining separated endpoints."""
    finite(tolerance)
    if tolerance <= 0 or not gaps:
        raise ValueError('BAD_DISTANCE_BAND_INPUT')
    for d in gaps.values():
        finite(d, nonnegative=True)
    ordered = sorted(gaps, key=lambda k: (gaps[k], k))
    bands = []
    for species in ordered:
        if not bands or gaps[species] - gaps[bands[-1][0]] > tolerance:
            bands.append([species])
        else:
            bands[-1].append(species)
    return bands


def band_ranks(gaps, tolerance):
    result = {}
    used = 0
    for block in distance_bands(gaps, tolerance):
        rank = used + (len(block)+1)/2
        result.update({s: rank for s in block})
        used += len(block)
    return result


def tier_ids(gaps, tolerance):
    tiers = {'near': [], 'middle': [], 'far': []}
    labels = tuple(tiers)
    used = 0
    for block in distance_bands(gaps, tolerance):
        midpoint = used + (len(block)-1)/2
        index = min(2, int(3*midpoint/len(gaps)))
        tiers[labels[index]].extend(block)
        used += len(block)
    return tiers


def midrank(values):
    for value in values:
        finite(value)
    order = sorted(range(len(values)), key=lambda i: (values[i], i))
    out = [0.0]*len(values)
    start = 0
    while start < len(order):
        end = start+1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        for i in range(start, end):
            out[order[i]] = (start+1+end)/2
        start = end
    return out


def pearson(xs, ys):
    if len(xs) != len(ys):
        raise ValueError('CORRELATION_LENGTH_MISMATCH')
    if len(xs) < 2:
        return None
    for val in list(xs)+list(ys):
        finite(val)
    mx = math.fsum(xs)/len(xs)
    my = math.fsum(ys)/len(ys)
    xx = math.fsum((v-mx)**2 for v in xs)
    yy = math.fsum((v-my)**2 for v in ys)
    if xx == 0 or yy == 0:
        return None
    return math.fsum((a-mx)*(b-my) for a,b in zip(xs,ys))/math.sqrt(xx*yy)


def spearman(xs, ys):
    return pearson(midrank(xs), midrank(ys))


def rnd(value, digits=8):
    return None if value is None else round(value, digits)


def neighbor_prediction(train, distances, endpoint, tolerance):
    if not train or set(train) != set(distances):
        raise ValueError('NEIGHBOR_TRAINING_ID_MISMATCH')
    finite(tolerance, nonnegative=True)
    if tolerance == 0:
        raise ValueError('ZERO_TOLERANCE_NOT_REGISTERED')
    for value in distances.values():
        finite(value, nonnegative=True)
    lowest = min(distances.values())
    tied = sorted(k for k in train if distances[k]-lowest <= tolerance)
    values = [finite(train[k]['responses'][endpoint]) for k in tied]
    return math.fsum(values)/len(values), tied


def expression_novelty(train, test_genes):
    if not train:
        raise ValueError('EMPTY_TRAINING_SET')
    xs = [[finite(row['genes'][g]) for g in BRANCH_GENES] for row in train.values()]
    mu = [math.fsum(x[j] for x in xs)/len(xs) for j in range(3)]
    sd = [math.sqrt(math.fsum((x[j]-mu[j])**2 for x in xs)/len(xs)) for j in range(3)]
    scale = [x if x > 1e-10 else 1.0 for x in sd]
    target = [finite(test_genes[g]) for g in BRANCH_GENES]
    return min(math.sqrt(math.fsum(((target[j]-x[j])/scale[j])**2 for j in range(3)))
               for x in xs)


def loo_rows(samples, pairs, models, predictor):
    """Pass only test genes (never test y) to the frozen predictor."""
    ids = sorted(samples)
    if len(ids) < 11 or not set(ids).issubset(pairs):
        raise ValueError('LOO_PANEL_OR_DISTANCE_COVERAGE')
    for a in ids:
        if pairs[a].get(a) != 0 or not set(ids).issubset(pairs[a]):
            raise ValueError('DISTANCE_MATRIX_DIAGONAL_OR_IDS')
        for b in ids:
            finite(pairs[a][b], nonnegative=True)
            if not math.isclose(pairs[a][b], pairs[b][a], rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError('DISTANCE_MATRIX_ASYMMETRY')
    rows = []
    for species in ids:
        train = {k:samples[k] for k in ids if k != species}
        test = {'genes': dict(samples[species]['genes'])}
        d = {k:pairs[species][k] for k in train}
        forecasts = {}
        neighbor_counts = {}
        for endpoint in ENDPOINTS:
            forecasts[endpoint] = {
                model: finite(predictor(train, test, endpoint, genes))
                for model,genes in models.items()
            }
            for tol in TOLERANCES:
                value, ties = neighbor_prediction(train,d,endpoint,tol)
                forecasts[endpoint]['phylo_nearest_'+tk(tol)] = value
                neighbor_counts[tk(tol)] = len(ties)
        rows.append({'id':species, 'n_train':len(train), 'gap':min(d.values()),
                     'novelty':expression_novelty(train,test['genes']),
                     'neighbor_counts':neighbor_counts, 'pred':forecasts,
                     'y':{e:finite(samples[species]['responses'][e]) for e in ENDPOINTS}})
    return rows


def losses(rows, endpoint, model):
    return [(r['y'][endpoint]-r['pred'][endpoint][model])**2 for r in rows]


def raw_rmse(rows, endpoint, model):
    if not rows:
        return None
    return math.sqrt(math.fsum(losses(rows,endpoint,model))/len(rows))


def risk(rows, endpoint, model):
    residuals = [r['y'][endpoint]-r['pred'][endpoint][model] for r in rows]
    if not residuals:
        return None
    return {'rmse':rnd(math.sqrt(math.fsum(x*x for x in residuals)/len(rows))),
            'mae':rnd(math.fsum(abs(x) for x in residuals)/len(rows))}


def gain(rows, endpoint, candidate, baseline):
    base = raw_rmse(rows,endpoint,baseline)
    return 100*(1-raw_rmse(rows,endpoint,candidate)/base) if base and base > 0 else None


def tol_summary(rows, endpoint, tolerance):
    gaps = {r['id']:r['gap'] for r in rows}
    ranks = band_ranks(gaps,tolerance)
    ranks_list = [ranks[r['id']] for r in rows]
    branch = losses(rows,endpoint,'branch_three')
    null = losses(rows,endpoint,'train_mean_baseline')
    difference = [a-b for a,b in zip(null,branch)]
    rho = spearman(ranks_list,branch)
    nearest = 'phylo_nearest_'+tk(tolerance)
    nearest_loss = losses(rows,endpoint,nearest)
    ids = {r['id']:r for r in rows}
    band_dist = []; band_error = []
    for group in distance_bands(gaps,tolerance):
        band_dist.append(math.fsum(gaps[s] for s in group)/len(group))
        band_error.append(math.fsum((ids[s]['y'][endpoint]-ids[s]['pred'][endpoint]['branch_three'])**2
                                    for s in group)/len(group))
    tiers = {}
    for label, group in tier_ids(gaps,tolerance).items():
        rr = [ids[s] for s in group]
        tiers[label] = None if not rr else {
            'n_species':len(rr), 'gap_min':rnd(min(gaps[s] for s in group)),
            'gap_max':rnd(max(gaps[s] for s in group)),
            'train_mean':risk(rr,endpoint,'train_mean_baseline'),
            'branch_three':risk(rr,endpoint,'branch_three'),
            'phylo_nearest':risk(rr,endpoint,nearest),
            'branch_rmse_gain_vs_mean_pct':rnd(gain(rr,endpoint,'branch_three','train_mean_baseline')),
            'branch_rmse_gain_vs_nearest_pct':rnd(gain(rr,endpoint,'branch_three',nearest)),
        }
    near,far = tiers['near'],tiers['far']
    pattern = None if rho is None or near is None or far is None else (
        rho > 0 and far['branch_three']['rmse'] > near['branch_three']['rmse'])
    a=risk(rows,endpoint,'branch_three'); b=risk(rows,endpoint,nearest)
    return {
        'n_distance_bands':len(band_dist),
        'rho_gap_branch_squared_error':rnd(rho),
        'rho_gap_mean_squared_error':rnd(spearman(ranks_list,null)),
        'rho_gap_mean_minus_branch_squared_error':rnd(spearman(ranks_list,difference)),
        'rho_equal_band_mean_gap_vs_mean_branch_squared_error':rnd(spearman(band_dist,band_error)),
        'branch_rmse_gain_vs_nearest_pct':rnd(gain(rows,endpoint,'branch_three',nearest)),
        'branch_beats_nearest_in_rmse_and_mae':a['rmse']<b['rmse'] and a['mae']<b['mae'],
        'n_species_branch_squared_error_better_than_nearest':sum(x<y for x,y in zip(branch,nearest_loss)),
        'distance_pattern_support':pattern,
        'tiers':tiers,
    }


def influence(rows,endpoint):
    """Drop already-scored rows, without any training refit or CI claim."""
    rhos=[]; gains=[]
    for i in range(len(rows)):
        subset = rows[:i]+rows[i+1:]
        gaps = {r['id']:r['gap'] for r in subset}
        ranks = band_ranks(gaps,PRIMARY_TOL)
        rho = spearman([ranks[r['id']] for r in subset],losses(subset,endpoint,'branch_three'))
        if rho is not None: rhos.append(rho)
        g = gain(subset,endpoint,'branch_three','phylo_nearest_'+tk(PRIMARY_TOL))
        if g is not None: gains.append(g)
    return {
        'already_scored_row_omission_not_refit_or_CI':True,
        'rho_gap_error_range': [rnd(min(rhos)),rnd(max(rhos))] if rhos else None,
        'branch_vs_nearest_rmse_gain_pct_range': [rnd(min(gains)),rnd(max(gains))] if gains else None,
        'n_omissions':len(rows),
    }


def summarize(rows, models):
    if len({r['id'] for r in rows}) != len(rows):
        raise ValueError('DUPLICATE_SCORED_SPECIES')
    gaps={r['id']:r['gap'] for r in rows}
    all_models=list(models)+['phylo_nearest_'+tk(t) for t in TOLERANCES]
    out={
        'n_species_scored_once':len(rows),
        'n_training_species_per_fold':len(rows)-1,
        'gap_min_median_max':[rnd(min(gaps.values())),rnd(median(gaps.values())),rnd(max(gaps.values()))],
        'distance_bands_by_tolerance':{
            tk(t):{'n_bands':len(distance_bands(gaps,t)),
                   'band_sizes':list(map(len,distance_bands(gaps,t))),
                   'nearest_training_tie_count_range':[min(r['neighbor_counts'][tk(t)] for r in rows),
                                                     max(r['neighbor_counts'][tk(t)] for r in rows)]}
            for t in TOLERANCES},
        'endpoints':{},
    }
    for endpoint in ENDPOINTS:
        branch=losses(rows,endpoint,'branch_three'); null=losses(rows,endpoint,'train_mean_baseline')
        out['endpoints'][endpoint]={
            'models':{name:risk(rows,endpoint,name) for name in all_models},
            'branch_rmse_gain_vs_mean_pct':rnd(gain(rows,endpoint,'branch_three','train_mean_baseline')),
            'n_species_branch_squared_error_better_than_mean':sum(a<b for a,b in zip(branch,null)),
            'raw_distance_rho_branch_squared_error':rnd(spearman([r['gap'] for r in rows],branch)),
            'rho_expression_novelty_branch_squared_error':rnd(spearman([r['novelty'] for r in rows],branch)),
            'rho_expression_novelty_mean_minus_branch_squared_error':rnd(spearman(
                [r['novelty'] for r in rows],[a-b for a,b in zip(null,branch)])),
            'tolerance_sensitivity':{tk(t):tol_summary(rows,endpoint,t) for t in TOLERANCES},
            'scored_row_influence':influence(rows,endpoint),
        }
    return out


def report(root):
    # These are existing parent modules, unchanged by v31.
    import analyze_original_petunieae_patristic_prediction_v30 as old30
    import analyze_original_petunieae_alternative_trees_v28 as old28
    import analyze_original_petunieae_2023_species_prediction_v26 as old26
    from analyze_original_petunieae_phylogeny_v27 import parse_newick
    from audit_petunieae_2023_processed_files_v26 import FILES, get_official_file, registered_location, verify_original
    paths={n:root/'results'/p for n,p in RECEIPTS.items()}
    old={n:json.loads(p.read_text()) for n,p in paths.items()}
    cache_original={};cache_published={}
    def source(url):
        if url not in cache_original: cache_original[url]=get_official_file(url)
        return cache_original[url]
    def published(url):
        if url not in cache_published: cache_published[url]=old28.fetch_published(url)
        return cache_published[url]
    replay=old30.analyze(**{'expected_v'+n:paths[n] for n in ('26','27','28','29')},
                         fetch_original=source,fetch_published=published)
    if replay != old['30']:
        raise ValueError('V30_ORIGINAL_SOURCE_OR_FROZEN_RESULT_DRIFT')
    samples,origin=old26.source_real_data(fetch=source)
    raw=source(registered_location(FILES['species_tree_pruned']))
    verify_original(raw,FILES['species_tree_pruned'])
    tree=parse_newick(raw.decode('utf-8-sig'))
    pairs=tree_distances(tree)
    reference=old30.patristic_pairs(tree)
    if set(pairs)!=set(samples)|{'BROW'} or len(samples)!=59:
        raise ValueError('ORIGINAL_60_TIP_59_INGROUP_GATE')
    for a in pairs:
        for b in pairs:
            if not math.isclose(pairs[a][b],reference[a][b],rel_tol=1e-12,abs_tol=1e-12):
                raise ValueError('INDEPENDENT_PATH_RECONSTRUCTION_MISMATCH')
    if tuple(old26.MODELS['branch_three']) != BRANCH_GENES or old26.RIDGE_LAMBDA != 10:
        raise ValueError('FIXED_MODEL_CHANGED')
    shared=old['29']['source_taxa_exact_codes']
    if len(shared)!=31 or len(set(shared))!=31 or not set(shared).issubset(samples):
        raise ValueError('FROZEN_31_FOCAL_ID_DRIFT')
    panels={}
    for label,ids in (('all_59_primary',sorted(samples)),('shared_31_secondary',shared)):
        subset={k:samples[k] for k in sorted(ids)}
        rows=loo_rows(subset,pairs,old26.MODELS,old26.ridge_predict)
        result=summarize(rows,old26.MODELS)
        if len(ids)==59:
            for e in ENDPOINTS:
                for model in old26.MODELS:
                    previous=old['26']['validation']['species_LOO']['endpoints'][e]['models'][model]
                    current=result['endpoints'][e]['models'][model]
                    if current['rmse']!=previous['out_of_fold_RMSE_log1p_mg_per_g'] or current['mae']!=previous['out_of_fold_MAE_log1p_mg_per_g']:
                        raise ValueError('V26_FROZEN_59_SPECIES_LOO_RMSE_MAE_DRIFT')
        panels[label]=result
    flags=[]
    for label,panel in panels.items():
        for e,obj in panel['endpoints'].items():
            primary=obj['tolerance_sensitivity'][tk(PRIMARY_TOL)]
            flags.append({'panel':label,'endpoint':e,'distance_pattern':primary['distance_pattern_support'],
                          'branch_beats_nearest_rmse_and_mae':primary['branch_beats_nearest_in_rmse_and_mae']})
    return {
        'schema':'anto_v31_source_verified_species_loo_distance_nearest_neighbor',
        'design':'EXPLORATORY_POST_V30_NOT_INDEPENDENT_BIOLOGICAL_REPLICATION',
        'parent_commit':PARENT,
        'original_chemical_study_doi':'10.1098/rspb.2023.0275',
        'original_source_tree_sha256':FILES['species_tree_pruned']['sha256'],
        'source_authentication_v26_to_v30_exact_replay':True,
        'frozen_v26_59_species_LOO_RMSE_and_MAE_reproduced_all_six_models':True,
        'decimal_paths_independently_agree_with_v30':True,
        'source_file_hashes':{k:{'md5':v['md5'],'sha256':v['sha256'],'bytes':v['size']}
                              for k,v in FILES.items() if k!='pigment_cluster_source'},
        'upstream_receipt_sha256':{n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in paths.items()},
        'primary_tolerance':PRIMARY_TOL,
        'mandatory_tolerance_sensitivity':list(TOLERANCES),
        'same_distance_bands_are_not_independent_lineages':True,
        'distance_units':'original_Newick_not_validated_years',
        'fixed_ridge_lambda':10,
        'panels':panels,
        'directional_diagnostics':flags,
        'all_four_distance_patterns_supported':all(x['distance_pattern'] is True for x in flags),
        'no_causal_p_values_or_new_functional_ancestral_regain_measurement':True,
        'raw_original_observations_not_redistributed':True,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--expected',type=Path)
    args=ap.parse_args()
    result=report(Path(__file__).resolve().parents[1])
    if args.expected and json.loads(args.expected.read_text())!=result:
        raise ValueError('V31_FROZEN_DERIVED_RESULT_DRIFT')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    brief={}
    for label,panel in result['panels'].items():
        brief[label]={'n':panel['n_species_scored_once'],
                      'distance_bands':panel['distance_bands_by_tolerance'], 'endpoints':{}}
        for e,obj in panel['endpoints'].items():
            p=obj['tolerance_sensitivity'][tk(PRIMARY_TOL)]
            brief[label]['endpoints'][e]={
                'models':obj['models'], 'branch_gain_vs_mean':obj['branch_rmse_gain_vs_mean_pct'],
                'branch_gain_vs_nearest':p['branch_rmse_gain_vs_nearest_pct'],
                'rho_gap_error':p['rho_gap_branch_squared_error'],
                'rho_gap_relative_gain':p['rho_gap_mean_minus_branch_squared_error'],
                'rho_molecular_novelty_error':obj['rho_expression_novelty_branch_squared_error'],
                'tiers':p['tiers'],'influence':obj['scored_row_influence']}
    print('ANTO_V31_SUMMARY',json.dumps({'panels':brief,'directional_diagnostics':result['directional_diagnostics']},sort_keys=True))
    full=json.dumps(result,sort_keys=True,separators=(',',':'))
    count=(len(full)+5999)//6000
    for i in range(count): print('ANTO_V31_CHUNK',i+1,count,full[6000*i:6000*(i+1)])

if __name__=='__main__': main()
