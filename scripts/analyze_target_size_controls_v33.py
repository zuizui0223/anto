"""v33: same-target, same-training-count controls for original clade holdouts.

Reuses v32 prediction primitives unchanged. New outcomes are conditional on the
same original 59-species study; random training subsets are not new studies.
"""
from __future__ import annotations
import argparse, hashlib, json, math, random
from collections import Counter
from pathlib import Path
import numpy as np

CAPS=(8,15)
REPS=31
SEED=20261114
MODELS=('mean','genes_only','phylo','phylo_genes')
ENDPOINTS=('anthocyanidins_log1p_total','flavonols_log1p_total')
PARENT='9185e9f2766df62df1a94b668f79841756263a38'


def control_train(ids,target,clade_size,cap,rep):
    """Outcome-blind training sample; target and sample count matched exactly."""
    ids=sorted(ids)
    if (len(ids)!=len(set(ids)) or target not in ids or cap not in CAPS or
        not isinstance(clade_size,int) or isinstance(clade_size,bool) or
        not 1<=clade_size<=cap or len(ids)-clade_size<4 or
        not isinstance(rep,int) or not 0<=rep<REPS):
        raise ValueError('INVALID_MATCHED_DESIGN')
    key=f'{SEED}|{cap}|{rep}|{target}'
    seed=int.from_bytes(hashlib.sha256(key.encode()).digest()[:8],'big')
    pool=[s for s in ids if s!=target]
    return sorted(random.Random(seed).sample(pool,len(ids)-clade_size))


def predict_subset(k,y,genes,train,test):
    """Unused taxa remain unused; only y[train] enters nested fitting."""
    import analyze_phylo_expression_v32 as old
    k=np.asarray(k,float);y=np.asarray(y,float);genes=np.asarray(genes,float)
    if y.ndim!=2 or y.shape[1]!=2:raise ValueError('ENDPOINT_SHAPE')
    n=len(y);tr=list(train);te=list(test)
    if (len(tr)<4 or not te or len(set(tr))!=len(tr) or len(set(te))!=len(te) or
        set(tr)&set(te) or not set(tr+te).issubset(range(n))):
        raise ValueError('TRAIN_TEST_IDENTITY')
    if k.shape!=(n,n) or genes.shape!=(n,3):raise ValueError('INPUT_SHAPE')
    kt=k[np.ix_(tr,tr)];kc=k[np.ix_(te,tr)]
    x,xt=old.design(genes[tr],genes[te],len(tr),len(te))
    pred={'mean':np.tile(y[tr].mean(0),(len(te),1)),
          'genes_only':xt@np.linalg.solve(x.T@x+np.diag([0,10,10,10]),x.T@y[tr]),
          'phylo':np.empty((len(te),2)), 'phylo_genes':np.empty((len(te),2))}
    selected=[]
    for e in range(2):
        eta,_=old.select_eta(kt,y[tr,e]);selected.append(eta)
        pred['phylo'][:,e]=old.predict(kt,kc,y[tr,e],eta)
        pred['phylo_genes'][:,e]=old.predict(kt,kc,y[tr,e],eta,genes[tr],genes[te])
    if any(not np.isfinite(p).all() for p in pred.values()):raise ValueError('NONFINITE_PREDICTION')
    return {'eta':selected,'pred':pred}


def risk_from_losses(squared,absolute):
    a=np.asarray(squared,float);b=np.asarray(absolute,float)
    if a.size==0 or a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(a<0) or np.any(b<0):
        raise ValueError('INVALID_LOSS')
    return {'rmse':round(float(np.sqrt(a.mean())),8),'mae':round(float(b.mean()),8)}


def paired_effect(real,controls,folds):
    """Difference of added-expression MSE gains, paired on original target."""
    rp=np.asarray(real['phylo'],float);rg=np.asarray(real['phylo_genes'],float)
    cp=np.asarray(controls['phylo'],float);cg=np.asarray(controls['phylo_genes'],float)
    if (rp.ndim!=1 or rg.shape!=rp.shape or cp.ndim!=2 or cg.shape!=cp.shape or
        cp.shape[1]!=len(rp) or cp.shape[0]<1 or
        sorted(i for f in folds for i in f)!=list(range(len(rp))) or any(not f for f in folds)):
        raise ValueError('PAIRED_LOSS_OR_FOLD_GEOMETRY')
    if any(not np.isfinite(a).all() or np.any(a<0) for a in (rp,rg,cp,cg)):
        raise ValueError('INVALID_PAIRED_LOSS')
    rgain=rp-rg;cgain=cp-cg
    delta=(rgain[None,:]-cgain).mean(0)  # subtract before averaging: identical controls stay exact zero
    ppen=float((rp-cp.mean(0)).mean());gpen=float((rg-cg.mean(0)).mean())
    if not math.isclose(float(delta.mean()),ppen-gpen,rel_tol=1e-10,abs_tol=1e-10):
        raise ValueError('LOSS_DECOMPOSITION_DRIFT')
    block=[float(delta[f].mean()) for f in folds]
    omission=[]
    for f in folds:
        keep=sorted(set(range(len(rp)))-set(f))
        if keep:omission.append(float(delta[keep].mean()))
    dist=(rgain.mean()-cgain.mean(1))
    return {'mean_increment_difference':round(float(delta.mean()),10),
            'equal_original_clade_increment_difference':round(float(np.mean(block)),10),
            'real_clade_mean_increment':round(float(rgain.mean()),10),
            'matched_control_mean_increment':round(float(cgain.mean()),10),
            'phylo_loss_penalty':round(ppen,10),'phylo_genes_loss_penalty':round(gpen,10),
            'n_positive_original_clades':sum(v>0 for v in block),'n_original_clades':len(folds),
            'all_original_clade_differences':[round(v,10) for v in block],
            'n_species_positive_difference':int(np.sum(delta>0)),
            'n_control_repetitions':cp.shape[0],
            'all_repetition_mean_differences':[round(float(v),10) for v in dist],
            'omit_one_scored_clade_mean_difference_range':
                [round(min(omission),10),round(max(omission),10)] if omission else None,
            'omission_is_not_refitting_or_CI':True}


def directional_gate(cells):
    return len(cells)==4 and all(c['mean_increment_difference']>0 and
        2*c['n_positive_original_clades']>c['n_original_clades'] for c in cells)


def evaluate(samples,ids,k,distances,fold_names,cap,previous):
    import analyze_phylo_expression_v32 as old
    index={s:i for i,s in enumerate(ids)}
    folds=[[index[s] for s in f] for f in fold_names]
    if len(ids)!=59 or sorted(i for f in folds for i in f)!=list(range(59)):
        raise ValueError('ORIGINAL_59_SPECIES_FOLD_GATE')
    y=np.array([[samples[s]['responses'][e] for e in ENDPOINTS] for s in ids])
    genes=np.array([[samples[s]['genes'][g] for g in old.GENES] for s in ids])
    actual={m:np.full((59,2),np.nan) for m in MODELS}
    controls={m:np.full((REPS,59,2),np.nan) for m in MODELS}
    group_for={i:f for f in folds for i in f}
    realgap=np.empty(59);randomgap=np.empty((REPS,59));peer_counts=np.empty((REPS,59),int)
    train_sizes=np.empty(59,int);eta_counts=[Counter(),Counter()];cache={};unique_controls=set()
    for fold in folds:
        tr=sorted(set(range(59))-set(fold))
        pred=predict_subset(k,y,genes,tr,fold)
        for m in MODELS:actual[m][fold]=pred['pred'][m]
        for j,i in enumerate(fold):
            realgap[i]=min(distances[ids[i]][ids[t]] for t in tr)
            train_sizes[i]=len(tr)
            cache[(tuple(tr),i)]={m:pred['pred'][m][j] for m in MODELS}
            cache[('eta',tuple(tr),i)]=pred['eta']
    for e,endpoint in enumerate(ENDPOINTS):
        for m in MODELS:
            got=old.metrics(y[:,e],actual[m][:,e]);wanted=previous['endpoints'][endpoint]['risks'][m]
            if any(not math.isclose(got[t],wanted[t],abs_tol=1e-8,rel_tol=0) for t in ('rmse','mae')):
                raise ValueError('V32_REAL_CLADE_MODEL_METRIC_DRIFT')
    for rep in range(REPS):
        for i,species in enumerate(ids):
            fold=group_for[i]
            names=control_train(ids,species,len(fold),cap,rep)
            tr=[index[s] for s in names]
            if len(tr)!=train_sizes[i] or i in tr:raise ValueError('MATCHING_INVARIANT_FAILURE')
            key=(tuple(tr),i);unique_controls.add(key)
            if key not in cache:
                p=predict_subset(k,y,genes,tr,[i])
                cache[key]={m:p['pred'][m][0] for m in MODELS}
                cache[('eta',tuple(tr),i)]=p['eta']
            for m in MODELS:controls[m][rep,i]=cache[key][m]
            for e,eta in enumerate(cache[('eta',tuple(tr),i)]):eta_counts[e][f'{eta:g}']+=1
            randomgap[rep,i]=min(distances[species][ids[t]] for t in tr)
            peer_counts[rep,i]=len(set(tr)&(set(fold)-{i}))
        print(f'V33_CONTROL_PROGRESS cap={cap} repetition={rep+1}/{REPS}',flush=True)
    if any(not np.isfinite(a).all() for a in list(actual.values())+list(controls.values())):
        raise ValueError('MISSING_FINITE_PREDICTIONS')
    singleton=[i for i in range(59) if len(group_for[i])==1]
    for i in singleton:
        for m in MODELS:
            if not np.array_equal(controls[m][:,i],np.tile(actual[m][i],(REPS,1))):
                raise ValueError('SINGLETON_MATCHED_CONTROL_NOT_IDENTICAL')
    out={'n_original_species':59,'n_repetitions':REPS,'n_original_clades':len(folds),
         'original_clade_source_tip_codes':fold_names,'train_count_histogram':dict(sorted(Counter(map(str,train_sizes)).items())),
         'same_target_and_training_count_verified':True,'n_singleton_targets':len(singleton),
         'n_unique_control_train_target_pairs':len(unique_controls),
         'distance_units':'original_Newick_not_validated_years',
         'training_composition_check':{
             'real_nearest_training_gap_mean':round(float(realgap.mean()),8),
             'control_nearest_training_gap_mean':round(float(randomgap.mean()),8),
             'mean_gap_difference_real_minus_control':round(float((realgap-randomgap).mean()),8),
             'fraction_control_gap_smaller_than_real':round(float((randomgap<realgap-1e-5).mean()),8),
             'mean_real_clade_peers_retained_in_control':round(float(peer_counts.mean()),8),
             'fraction_controls_with_any_real_clade_peer':round(float((peer_counts>0).mean()),8)},
         'endpoints':{}}
    for e,endpoint in enumerate(ENDPOINTS):
        rsq={m:(y[:,e]-actual[m][:,e])**2 for m in MODELS}
        csq={m:(y[:,e]-controls[m][:,:,e])**2 for m in MODELS}
        real_risk={m:risk_from_losses(rsq[m],np.abs(y[:,e]-actual[m][:,e])) for m in MODELS}
        control_risk={m:risk_from_losses(csq[m],np.abs(y[:,e]-controls[m][:,:,e])) for m in MODELS}
        repeated=[{m:old.metrics(y[:,e],controls[m][rep,:,e]) for m in MODELS} for rep in range(REPS)]
        def increment(risk):return round(100*(1-risk['phylo_genes']['rmse']/risk['phylo']['rmse']),6)
        out['endpoints'][endpoint]={
            'real_clade_risks':real_risk,'matched_control_pooled_risks':control_risk,
            'real_clade_rmse_gain_pct':increment(real_risk),
            'matched_control_pooled_rmse_gain_pct':increment(control_risk),
            'all_control_repetition_risks':repeated,
            'control_selected_eta_counts':dict(sorted(eta_counts[e].items())),
            'paired_effect':paired_effect(rsq,csq,folds)}
    return out


def report(root,cache_dir):
    import analyze_phylo_expression_v32 as old
    receipt=json.loads((root/'results/petunieae_phylo_expression_v32_20261114.receipt.json').read_text())
    replay=old.report(root,cache_dir)
    serialized=(json.dumps(replay,sort_keys=True,indent=2)+'\n').encode()
    if hashlib.sha256(serialized).hexdigest()!=receipt['full_result_sha256']:
        raise ValueError('V32_FULL_SOURCE_REPLAY_DIGEST_MISMATCH')
    import audit_petunieae_2023_processed_files_v26 as source
    from analyze_original_petunieae_2023_species_prediction_v26 import source_real_data
    from analyze_original_petunieae_phylogeny_v27 import parse_newick,prune_outgroup,exact_monophyletic_partitions
    from analyze_petunieae_loo_distance_v31 import tree_distances
    registered={source.registered_location(v):v for key,v in source.FILES.items() if key!='pigment_cluster_source'}
    def verified(url):
        record=registered[url];raw=(cache_dir/(record['file_id']+'.raw')).read_bytes()
        source.verify_original(raw,record);return raw
    samples,_=source_real_data(fetch=verified)
    raw=verified(source.registered_location(source.FILES['species_tree_pruned']))
    tree=parse_newick(raw.decode('utf-8-sig'));tip_ids,kfull=old.shared_covariance(tree)
    if set(tip_ids)!=set(samples)|{'BROW'} or len(samples)!=59:raise ValueError('SOURCE_SPECIES_MISMATCH')
    ids=sorted(samples);tipix={s:i for i,s in enumerate(tip_ids)};idx=[tipix[s] for s in ids]
    k=kfull[np.ix_(idx,idx)];dist=tree_distances(tree)
    for i,a in enumerate(ids):
        for j,b in enumerate(ids):
            if not math.isclose(k[i,i]+k[j,j]-2*k[i,j],dist[a][b],abs_tol=1e-10,rel_tol=1e-10):
                raise ValueError('COVARIANCE_DISTANCE_MISMATCH')
    designs={}
    for cap in CAPS:
        key=f'59_species_clade{cap}'
        folds=exact_monophyletic_partitions(prune_outgroup(tree),cap)
        designs[key]=evaluate(samples,ids,k,dist,folds,cap,replay['designs'][key])
    cells=[obj['paired_effect'] for design in designs.values() for obj in design['endpoints'].values()]
    return {'schema':'anto_v33_same_target_training_count_composition_control',
        'parent_commit':PARENT,'source_DOI':'10.1098/rspb.2023.0275',
        'v32_full_original_source_replayed_sha256':receipt['full_result_sha256'],
        'source_tree_sha256':source.FILES['species_tree_pruned']['sha256'],
        'matched_target_n':59,'control_repetitions':REPS,'seed_base':SEED,
        'unchanged_models':{'nuggets':list(old.ETAS),'gene_lambda':old.GENE_LAMBDA,'genes':list(old.GENES)},
        'designs':designs,'all_four_directional_composition_checks':directional_gate(cells),
        'primary_estimand':'clade molecular MSE increment minus same-target/count control mean increment',
        'controls_not_independent_studies_or_permutation_p_values':True,
        'random_training_changes_more_than_phylogenetic_distance':True,
        'no_new_original_chemical_study_or_causal_ancestral_regain':True,
        'v32_primary_failure_not_overridden':True}


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-cache',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--expected',type=Path)
    args=p.parse_args();r=report(Path(__file__).resolve().parents[1],args.source_cache)
    if args.expected and json.loads(args.expected.read_text())!=r:raise ValueError('V33_FROZEN_RESULT_DRIFT')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    short=[]
    for name,design in r['designs'].items():
        for ep,v in design['endpoints'].items():
            short.append({'design':name,'endpoint':ep,'real_risks':v['real_clade_risks'],
                'control_risks':v['matched_control_pooled_risks'],
                'real_rmse_gain_pct':v['real_clade_rmse_gain_pct'],
                'control_rmse_gain_pct':v['matched_control_pooled_rmse_gain_pct'],
                'paired_effect':v['paired_effect'],'composition_check':design['training_composition_check']})
    print('ANTO_V33_SUMMARY',json.dumps({'comparisons':short,'gate':r['all_four_directional_composition_checks']},sort_keys=True))
    print('ANTO_V33_RESULT_SHA256',hashlib.sha256(args.output.read_bytes()).hexdigest())

if __name__=='__main__':main()
