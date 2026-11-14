"""Nested training-only phylogenetic prediction versus the same model + expression.

A single inspected 2023 chemistry dataset, not causal or independent replication.
Original sources and the entire v31 receipt must replay before new calculations.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor, cho_solve

ETAS=(.01,.1,1.,10.)
GENE_LAMBDA=10.
GENES=('F3primeH','F35H','FLS')
ENDPOINTS=('anthocyanidins_log1p_total','flavonols_log1p_total')
PARENT='45154b7fb4ad653e6c380afbc2a707df40c9a6fa'


def shared_covariance(root):
    paths={};seen=set()
    def visit(node,path):
        if id(node) in seen:raise ValueError('TREE_REUSED_NODE')
        seen.add(id(node))
        if node.children:
            for child in node.children:
                edge=child.length
                if edge is None or not math.isfinite(edge) or edge<0:raise ValueError('INVALID_SOURCE_EDGE')
                visit(child,path+[(id(child),float(edge))])
        else:
            if not node.label or node.label in paths:raise ValueError('BAD_SOURCE_TIP')
            paths[node.label]=dict(path)
    visit(root,[])
    ids=sorted(paths)
    if len(ids)<2:raise ValueError('TOO_FEW_SOURCE_TIPS')
    k=np.array([[math.fsum(paths[a][e] for e in paths[a].keys()&paths[b].keys())
                 for b in ids] for a in ids],dtype=float)
    if not np.isfinite(k).all() or not np.any(np.diag(k)>0):raise ValueError('INVALID_SHARED_PATH_COVARIANCE')
    if np.linalg.eigvalsh(k).min() < -1e-8:raise ValueError('SOURCE_COVARIANCE_NOT_PSD')
    return ids,k


def design(gtrain,gtest,n,ntest):
    if gtrain is None and gtest is None:return np.ones((n,1)),np.ones((ntest,1))
    a=np.asarray(gtrain,float);b=np.asarray(gtest,float)
    if a.shape!=(n,3) or b.shape!=(ntest,3) or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('GENE_COORDINATES_INVALID')
    mu=a.mean(0);sd=a.std(0);sd=np.where(sd>1e-10,sd,1.)
    return np.column_stack((np.ones(n),(a-mu)/sd)),np.column_stack((np.ones(ntest),(b-mu)/sd))


def predict(ktrain,kcross,ytrain,eta,gtrain=None,gtest=None):
    """Penalized GLS trend + conditional phylogenetic residual mean."""
    k=np.asarray(ktrain,float);c=np.asarray(kcross,float);y=np.asarray(ytrain,float)
    n=len(y)
    if y.ndim!=1 or n<2 or k.shape!=(n,n) or c.ndim!=2 or c.shape[1]!=n:
        raise ValueError('PREDICTION_SHAPE')
    if not np.isfinite(k).all() or not np.isfinite(c).all() or not np.isfinite(y).all():
        raise ValueError('NONFINITE_TRAINING')
    if not np.allclose(k,k.T,rtol=1e-12,atol=1e-12):raise ValueError('ASYMMETRIC_COVARIANCE')
    if eta not in ETAS:raise ValueError('UNDECLARED_NUGGET')
    scale=np.diag(k).mean()
    if not scale>0:raise ValueError('ZERO_COVARIANCE_SCALE')
    covariance=k/scale+eta*np.eye(n)
    factor=cho_factor(covariance,lower=True,check_finite=True)
    x,xt=design(gtrain,gtest,n,len(c))
    invx=cho_solve(factor,x);invy=cho_solve(factor,y)
    penalty=np.diag([0.]+[GENE_LAMBDA]*(x.shape[1]-1))
    beta=np.linalg.solve(x.T@invx+penalty,x.T@invy)
    answer=xt@beta+(c/scale)@cho_solve(factor,y-x@beta)
    if not np.isfinite(answer).all():raise ValueError('NONFINITE_PREDICTION')
    return answer


def select_eta(ktrain,ytrain):
    """Inner species LOO of P ONLY; no outer test outcomes/features."""
    k=np.asarray(ktrain,float);y=np.asarray(ytrain,float);n=len(y)
    if n<4 or k.shape!=(n,n):raise ValueError('INNER_TRAINING_TOO_SMALL')
    scores=[]
    for eta in ETAS:
        residual=[]
        for i in range(n):
            keep=[j for j in range(n) if j!=i]
            value=predict(k[np.ix_(keep,keep)],k[np.ix_([i],keep)],y[keep],eta)[0]
            residual.append((y[i]-value)**2)
        scores.append(float(np.mean(residual)))
    return ETAS[int(np.argmin(scores))],scores


def outer_fold(k,y,genes,train,test):
    k=np.asarray(k,float);y=np.asarray(y,float);genes=np.asarray(genes,float)
    if y.ndim!=2:raise ValueError('ENDPOINT_MATRIX_REQUIRED')
    n=len(y);tr=list(train);te=list(test)
    if len(set(tr))!=len(tr) or len(set(te))!=len(te) or set(tr)&set(te) or set(tr)|set(te)!=set(range(n)):
        raise ValueError('OUTER_FOLD_IDENTITY_OR_LEAKAGE')
    if k.shape!=(n,n) or genes.shape!=(n,3) or len(tr)<4 or not te:raise ValueError('OUTER_SHAPE')
    kt=k[np.ix_(tr,tr)];kc=k[np.ix_(te,tr)]
    x,xt=design(genes[tr],genes[te],len(tr),len(te))
    pred={'mean':np.tile(y[tr].mean(0),(len(te),1)),
          'genes_only':xt@np.linalg.solve(x.T@x+np.diag([0,10,10,10]),x.T@y[tr]),
          'phylo':np.zeros((len(te),y.shape[1])),
          'phylo_genes':np.zeros((len(te),y.shape[1]))}
    for eta in ETAS:
        pred[f'phylo_eta{eta:g}']=np.zeros((len(te),y.shape[1]))
        pred[f'phylo_genes_eta{eta:g}']=np.zeros((len(te),y.shape[1]))
    selected=[];inner=[]
    for e in range(y.shape[1]):
        eta,scores=select_eta(kt,y[tr,e]);selected.append(eta);inner.append(scores)
        for fixed in ETAS:
            a=predict(kt,kc,y[tr,e],fixed)
            b=predict(kt,kc,y[tr,e],fixed,genes[tr],genes[te])
            pred[f'phylo_eta{fixed:g}'][:,e]=a
            pred[f'phylo_genes_eta{fixed:g}'][:,e]=b
            if fixed==eta:pred['phylo'][:,e]=a;pred['phylo_genes'][:,e]=b
    return {'eta':selected,'inner_risks':inner,'pred':pred}


def metrics(y,p):
    residual=np.asarray(y)-np.asarray(p)
    return {'rmse':round(float(np.sqrt(np.mean(residual**2))),8),
            'mae':round(float(np.mean(np.abs(residual))),8)}


def improves(candidate,baseline):
    return candidate['rmse']<baseline['rmse'] and candidate['mae']<baseline['mae']


def contrast(y,base,candidate,folds):
    a=(y-base)**2;b=(y-candidate)**2;ra=float(np.sqrt(a.mean()));rb=float(np.sqrt(b.mean()))
    omission=[]
    for i in range(len(y)):
        keep=np.arange(len(y))!=i;den=float(np.sqrt(a[keep].mean()))
        if den>0:omission.append(100*(1-float(np.sqrt(b[keep].mean()))/den))
    block_gains=[float((a[f]-b[f]).mean()) for f in folds]
    return {'RMSE_reduction_pct':round(100*(1-rb/ra),6) if ra else None,
            'mean_squared_error_gain':round(float((a-b).mean()),8),
            'n_species_lower_squared_error':int(np.sum(b<a)),
            'n_species':len(y),'n_outer_folds':len(folds),
            'equal_fold_MSE_gain':round(float(np.mean(block_gains)),8),
            'n_folds_positive':sum(v>0 for v in block_gains),
            'all_fold_MSE_gains':[round(v,8) for v in block_gains],
            'omit_one_scored_species_gain_pct_range':[round(min(omission),6),round(max(omission),6)] if omission else None,
            'omission_is_not_refit_or_CI':True}


def evaluate(samples,ids,k,fold_names):
    ix={s:i for i,s in enumerate(ids)}
    folds=[[ix[s] for s in f] for f in fold_names]
    if sorted(i for f in folds for i in f)!=list(range(len(ids))):raise ValueError('SPECIES_NOT_SCORED_EXACTLY_ONCE')
    y=np.array([[samples[s]['responses'][e] for e in ENDPOINTS] for s in ids])
    x=np.array([[samples[s]['genes'][g] for g in GENES] for s in ids])
    predictions={};etas=[];inner_risks=[]
    for f in folds:
        tr=sorted(set(range(len(ids)))-set(f))
        result=outer_fold(k,y,x,tr,f);etas.append(result['eta']);inner_risks.append(result['inner_risks'])
        for name,vals in result['pred'].items():
            if name not in predictions:predictions[name]=np.full_like(y,np.nan)
            predictions[name][f]=vals
    if any(not np.isfinite(p).all() for p in predictions.values()):raise ValueError('MISSING_FINITE_OOF')
    out={'n_species':len(ids),'n_outer_folds':len(folds),'fold_sizes':[len(f) for f in folds],
         'endpoints':{}}
    for e,label in enumerate(ENDPOINTS):
        risks={name:metrics(y[:,e],p[:,e]) for name,p in predictions.items()}
        out['endpoints'][label]={
            'risks':risks,
            'selected_eta_counts':dict(sorted(Counter(f'{z[e]:g}' for z in etas).items())),
            'mean_inner_phylo_risks_by_eta':dict(zip(map(lambda z:f'{z:g}',ETAS),
                 [round(float(v),8) for v in np.mean([r[e] for r in inner_risks],axis=0)])),
            'added_expression':contrast(y[:,e],predictions['phylo'][:,e],predictions['phylo_genes'][:,e],folds),
            'added_expression_beats_phylo_RMSE_and_MAE':improves(risks['phylo_genes'],risks['phylo']),
            'fixed_eta_sensitivity':{f'{eta:g}':contrast(y[:,e],predictions[f'phylo_eta{eta:g}'][:,e],
                 predictions[f'phylo_genes_eta{eta:g}'][:,e],folds) for eta in ETAS}}
    return out


def report(root,cache_dir):
    import audit_petunieae_2023_processed_files_v26 as source
    registered={source.registered_location(v):v for key,v in source.FILES.items() if key!='pigment_cluster_source'}
    def verified(url):
        if url not in registered:raise ValueError('UNREGISTERED_SOURCE_REQUEST')
        record=registered[url];raw=(cache_dir/(record['file_id']+'.raw')).read_bytes()
        source.verify_original(raw,record);return raw
    source.get_official_file=verified
    import analyze_petunieae_loo_distance_v31 as v31
    original=root/'results/petunieae_species_loo_distance_v31_20261114.json'
    frozen=json.loads(original.read_text())
    if v31.report(root)!=frozen:raise ValueError('V31_SOURCE_REPLAY_MISMATCH')
    from analyze_original_petunieae_2023_species_prediction_v26 import source_real_data
    from analyze_original_petunieae_phylogeny_v27 import parse_newick,prune_outgroup,exact_monophyletic_partitions
    samples,_=source_real_data(fetch=verified)
    raw=verified(source.registered_location(source.FILES['species_tree_pruned']))
    tree=parse_newick(raw.decode('utf-8-sig'));tip_ids,fullk=shared_covariance(tree)
    if set(tip_ids)!=set(samples)|{'BROW'} or len(samples)!=59:raise ValueError('ORIGINAL_TAXON_GATE')
    distances=v31.tree_distances(tree)
    for i,a in enumerate(tip_ids):
        for j,b in enumerate(tip_ids):
            if not math.isclose(fullk[i,i]+fullk[j,j]-2*fullk[i,j],distances[a][b],abs_tol=1e-10,rel_tol=1e-10):
                raise ValueError('SHARED_PATH_VERSUS_PATRISTIC_DISAGREEMENT')
    old29=json.loads((root/'results/petunieae_fine_discordance_label_reference_v29_20261114.json').read_text())
    shared=old29['source_taxa_exact_codes']
    if len(set(shared))!=31 or not set(shared).issubset(samples):raise ValueError('SECONDARY_PANEL_GATE')
    allids=sorted(samples);subids=sorted(shared);ingroup=prune_outgroup(tree)
    definitions={
        '59_species_LOO_primary':(allids,[[s] for s in allids]),
        '31_species_LOO_secondary':(subids,[[s] for s in subids]),
        '59_species_clade8':(allids,exact_monophyletic_partitions(ingroup,8)),
        '59_species_clade15':(allids,exact_monophyletic_partitions(ingroup,15))}
    designs={};tipix={s:i for i,s in enumerate(tip_ids)}
    for name,(ids,folds) in definitions.items():
        idx=[tipix[s] for s in ids];designs[name]=evaluate(samples,ids,fullk[np.ix_(idx,idx)],folds)
        if name in ('59_species_LOO_primary','31_species_LOO_secondary'):
            prior=frozen['panels']['all_59_primary' if len(ids)==59 else 'shared_31_secondary']
            for endpoint in ENDPOINTS:
                for current,old in [('mean','train_mean_baseline'),('genes_only','branch_three')]:
                    for metric in ('rmse','mae'):
                        a=designs[name]['endpoints'][endpoint]['risks'][current][metric]
                        b=prior['endpoints'][endpoint]['models'][old][metric]
                        if not math.isclose(a,b,abs_tol=1e-8,rel_tol=0):raise ValueError('PARENT_LOO_COMPARATOR_DRIFT')
    return {'schema':'anto_v32_nested_phylogenetic_residual_prediction_gene_increment',
        'parent_commit':PARENT,'source_DOI':'10.1098/rspb.2023.0275',
        'v31_and_upstream_original_source_replay_exact':True,
        'v31_receipt_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),
        'source_tree_sha256':source.FILES['species_tree_pruned']['sha256'],
        'shared_path_covariance_independently_validated':True,
        'fixed_nuggets':list(ETAS),'fixed_gene_penalty':GENE_LAMBDA,
        'nugget_selection':'inner training species LOO of phylogeny-only; same selected eta in P and P+G',
        'designs':designs,
        'primary_both_endpoints_improve_rmse_and_mae':all(d['added_expression_beats_phylo_RMSE_and_MAE']
            for d in designs['59_species_LOO_primary']['endpoints'].values()),
        'exploratory_one_original_study_not_independent_biological_replication':True,
        'Brownian_covariance_is_an_assumption_not_verified_evolutionary_process':True,
        'no_causal_gene_function_or_ancestral_regain_or_descendant_fitness':True}


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-cache',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--expected',type=Path)
    args=p.parse_args();result=report(Path(__file__).resolve().parents[1],args.source_cache)
    if args.expected and json.loads(args.expected.read_text())!=result:raise ValueError('V32_FROZEN_OUTPUT_DRIFT')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    compact={name:{e:{'mean':d['risks']['mean'],'genes':d['risks']['genes_only'],
        'phylo':d['risks']['phylo'],'phylo_genes':d['risks']['phylo_genes'],
        'gain':d['added_expression'],'eta_counts':d['selected_eta_counts'],
        'fixed_eta_gain_pct':{k:v['RMSE_reduction_pct'] for k,v in d['fixed_eta_sensitivity'].items()}}
        for e,d in row['endpoints'].items()} for name,row in result['designs'].items()}
    print('ANTO_V32_SUMMARY',json.dumps(compact,sort_keys=True,separators=(',',':')))
    print('ANTO_V32_VERIFIED_FILE_SHA256',hashlib.sha256(args.output.read_bytes()).hexdigest())
    full=json.dumps(result,sort_keys=True,separators=(',',':'));count=(len(full)+5999)//6000
    for i in range(count):print('ANTO_V32_CHUNK',i+1,count,full[i*6000:(i+1)*6000])

if __name__=='__main__':main()
