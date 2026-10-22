import collections,itertools,json,math,shutil,subprocess
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from .common import RANKS,require,read_tsv,write_tsv,atomic_json,sha,digest,checked_checkpoint,commit,dependency_versions,r_versions
from .config import resolve

EXTRA={'Feature','Source','Trait','Status'}|set(RANKS)
def block_matrix(path,kind):
    rows=read_tsv(path);require(rows,f'Empty block: {path}');samples=[k for k in rows[0] if k not in EXTRA];features=[r['Feature'] for r in rows];require(len(features)==len(set(features)),'Duplicate feature IDs')
    X=np.array([[float(r[s]) if r[s] not in {'NA','','NaN'} else np.nan for r in rows] for s in samples]);require(np.all(X[np.isfinite(X)]>=0),'Negative abundances')
    if kind=='taxonomy':
        totals=X.sum(axis=1);X=np.divide(X,totals[:,None],out=np.full_like(X,np.nan),where=totals[:,None]>0)
    labels={r['Feature']:(r.get('Trait') or next((r[k] for k in reversed(RANKS[:-1]) if r.get(k) and not r[k].startswith('Unclassified')),r['Feature'])) for r in rows}
    return samples,features,X,labels

def matrix_tsv(path,D,names):write_tsv(path,['Sample']+names,[dict(Sample=s,**{t:float(D[i,j]) if np.isfinite(D[i,j]) else 'NA' for j,t in enumerate(names)}) for i,s in enumerate(names)])
def selected_meta(task,metadata):
    result=[]
    for r in metadata:
        keep=True
        for col,val in task.get('subset',{}).items():
            allowed=val if isinstance(val,list) else [val]
            if r.get(col) not in allowed:keep=False
        if keep:result.append(r)
    return result

def prepare_jobs(cfg,tables,functions):
    blocks=[(n,p,'taxonomy') for n,p in tables]+[(n,p,'traits') for n,p in functions];cache={n:(p,k,block_matrix(p,k) if read_tsv(p) else None) for n,p,k in blocks};metadata=read_tsv(resolve(cfg,cfg['metadata']));tasks=json.loads(resolve(cfg,cfg['tasks']).read_text());out=Path(cfg['_out'])/'statistics';out.mkdir(parents=True,exist_ok=True);jobs=[];summaries=[]
    for t in tasks:
        md=selected_meta(t,metadata);require(md,f'Empty task subset {t["id"]}')
        require(all(r.get(t['group'],'') for r in md),f'Missing group values in {t["id"]}')
        for block in t.get('blocks') or list(cache):
            require(block in cache,f'Task {t["id"]} requests unconfigured block {block}');path,kind,payload=cache[block];folder=out/t['id']/block;folder.mkdir(parents=True,exist_ok=True)
            if payload is None:
                job={'id':t['id']+'__'+block,'task':t['id'],'block':block,'kind':kind,'folder':str(folder),'metadata':str(folder/'samples.tsv'),'group':t['group'],'covariates':t.get('covariates',[]),'strata':t.get('strata'),'inferential_status':'not_estimable_empty_block','status_only':True,'missing_samples':[],'empty_samples':[],'N':len(md),'groups':dict(collections.Counter(r[t['group']] for r in md)),'note':'No target features or no positive functional reference annotations; not evidence of biological absence.','seed':cfg.get('seed',20261021),'permutations':cfg['statistics'].get('permutations',9999)}
                write_tsv(folder/'samples.tsv',list(md[0]),md);atomic_json(folder/'design.json',job);jobs.append(job);continue
            all_ids,featureids,fullX,labels=payload;m=[r for r in md if r['sample_id'] in all_ids];missing=[r['sample_id'] for r in md if r['sample_id'] not in all_ids]
            ids=[r['sample_id'] for r in m];require(ids,f'No available samples in {t["id"]}/{block}');X=fullX[[all_ids.index(s) for s in ids]];valid=np.all(np.isfinite(X),axis=1);empty=[ids[j] for j in np.flatnonzero(~valid)];m=[r for j,r in enumerate(m) if valid[j]];ids=[r['sample_id'] for r in m];X=X[valid];require(len(ids)>0,f'No non-empty samples in {t["id"]}/{block}')
            groups=[r[t['group']] for r in m];levels=sorted(set(groups));group_sizes=collections.Counter(groups);min_n=min(group_sizes.values());write_tsv(folder/'samples.tsv',list(m[0]),m)
            # Prevalence is evaluated globally per block, not separately per group.
            keep=np.any(np.nan_to_num(fullX)>0,axis=0)
            if kind=='taxonomy':keep=np.mean(np.nan_to_num(fullX)>=cfg['statistics'].get('relative_threshold',.0001),axis=0)>=cfg['statistics'].get('prevalence',.1)
            X=X[:,keep];fids=[v for v,k in zip(featureids,keep) if k];require(len(fids)>0,f'No retained features in {block}')
            write_tsv(folder/'feature_labels.tsv',['Feature','Label'],[dict(Feature=f,Label=labels.get(f,f)) for f in fids])
            write_tsv(folder/'matrix.tsv',['Sample']+fids,[dict(Sample=s,**dict(zip(fids,X[j]))) for j,s in enumerate(ids)])
            native=cdist(X,X,'braycurtis') if kind=='taxonomy' else cdist(np.sqrt(X),np.sqrt(X),'euclidean');bray=cdist(X,X,'braycurtis');np.fill_diagonal(bray,0);np.fill_diagonal(native,0);matrix_tsv(folder/'native_distance.tsv',native,ids);matrix_tsv(folder/'Bray_Curtis.tsv',bray,ids)
            summary=[];means=np.full((len(levels),len(levels)),np.nan)
            for i,a in enumerate(levels):
                ia=np.flatnonzero(np.array(groups)==a)
                for j,b in enumerate(levels):
                    ib=np.flatnonzero(np.array(groups)==b);vals=bray[np.ix_(ia,ib)];vals=vals[np.triu_indices(len(ia),1)] if a==b else vals.ravel();vals=vals[np.isfinite(vals)];means[i,j]=np.mean(vals) if len(vals) else np.nan
                    if j>=i:
                        record=dict(Task=t['id'],Block=block,Group_A=a,Group_B=b,N_A=len(ia),N_B=len(ib),Distance_Pairs=len(vals),Mean=float(np.mean(vals)) if len(vals) else 'NA',Median=float(np.median(vals)) if len(vals) else 'NA',Minimum=float(np.min(vals)) if len(vals) else 'NA',Maximum=float(np.max(vals)) if len(vals) else 'NA',SD_Descriptive=float(np.std(vals,ddof=1)) if len(vals)>1 else 'NA',Interpretation='Descriptive; pair distances are dependent');summary.append(record);summaries.append(record)
            write_tsv(folder/'all_group_pairs.tsv',list(summary[0]),summary);matrix_tsv(folder/'mean_pairwise_Bray_Curtis.tsv',means,levels)
            centers=np.stack([X[np.array(groups)==g].mean(axis=0) for g in levels]);centerD=cdist(centers,centers,'braycurtis');np.fill_diagonal(centerD,0);matrix_tsv(folder/'group_mean_profile_Bray_Curtis.tsv',centerD,levels)
            write_tsv(folder/'all_sample_pairs.tsv',['Sample_A','Sample_B','Group_A','Group_B','Bray_Curtis'],[dict(Sample_A=ids[i],Sample_B=ids[j],Group_A=groups[i],Group_B=groups[j],Bray_Curtis=bray[i,j] if np.isfinite(bray[i,j]) else 'NA') for i in range(len(ids)) for j in range(i+1,len(ids))])
            if kind=='taxonomy' and cfg.get('counts_provenance','observed')=='observed':
                original=read_tsv(path);alpha=[];rare=[];depth=cfg['statistics'].get('rarefaction_depth',1000);repeats=cfg['statistics'].get('rarefaction_repeats',200);rng=np.random.default_rng(cfg.get('seed',20261021))
                for sample,group in zip(ids,groups):
                    counts=np.array([float(r[sample]) for r in original]);require(np.all(counts==np.floor(counts)),'Alpha diversity requires integer observed counts, not fractions');counts=counts.astype(np.int64);total=int(counts.sum());positive=counts[counts>0];p=positive/total if total else np.array([]);f1=int(sum(counts==1));f2=int(sum(counts==2));alpha.append(dict(Sample=sample,Group=group,Reads=total,Observed_Feature_Rows=len(positive),Shannon=float(-np.sum(p*np.log(p))),Simpson=float(1-np.sum(p*p)) if total else 'NA',Chao1_Bias_Corrected=len(positive)+f1*(f1-1)/(2*(f2+1)),Goods_Coverage=1-f1/total if total else 'NA'))
                    if total>=depth:
                        values=[]
                        for _ in range(repeats):
                            z=rng.multivariate_hypergeometric(counts,depth);pp=z[z>0]/depth;values.append([np.sum(z>0),-np.sum(pp*np.log(pp)),1-np.sum(pp*pp)])
                        arr=np.array(values);rare.append(dict(Sample=sample,Group=group,Depth=depth,Repeats=repeats,Observed_Mean=float(arr[:,0].mean()),Shannon_Mean=float(arr[:,1].mean()),Simpson_Mean=float(arr[:,2].mean()),Status='Monte_Carlo_means_per_biological_sample'))
                    else:rare.append(dict(Sample=sample,Group=group,Depth=depth,Repeats=0,Observed_Mean='NA',Shannon_Mean='NA',Simpson_Mean='NA',Status='below_depth_not_scaled_or_duplicated'))
                write_tsv(folder/'alpha_observed_counts.tsv',list(alpha[0]),alpha);write_tsv(folder/'alpha_rarefied.tsv',list(rare[0]),rare)
            status='eligible_for_inference' if len(levels)>=2 and min_n>=2 and not missing and not empty and np.all(np.isfinite(native)) else 'descriptive_only'
            job={'id':t['id']+'__'+block,'task':t['id'],'block':block,'kind':kind,'folder':str(folder),'metadata':str(folder/'samples.tsv'),'group':t['group'],'covariates':t.get('covariates',[]),'strata':t.get('strata'),'inferential_status':status,'missing_samples':missing,'empty_samples':empty,'N':len(ids),'groups':dict(group_sizes),'note':t.get('note',''),'seed':cfg.get('seed',20261021),'permutations':cfg['statistics'].get('permutations',9999)};atomic_json(folder/'design.json',job);jobs.append(job)
    atomic_json(out/'jobs.json',jobs);write_tsv(out/'ALL_GROUP_PAIR_COMPARISONS.tsv',list(summaries[0]) if summaries else ['Task','Block'],summaries);return jobs

def run_statistics(cfg,tables,functions):
    root=Path(cfg['_out'])/'statistics'
    signature=digest({'algorithm':'ecology-v1','dependencies':dependency_versions(['numpy','scipy']),'R_runtime':r_versions(['vegan','permute','jsonlite']),'Python_code':sha(Path(__file__)),'seed':cfg.get('seed',20261021),'sources':{n:sha(p) for n,p in tables+functions},'metadata':sha(resolve(cfg,cfg['metadata'])),'tasks':sha(resolve(cfg,cfg['tasks'])),'statistics':cfg['statistics'],'R_code':sha(Path(__file__).with_name('statistics.R'))})
    if checked_checkpoint(root,signature):return json.loads((root/'jobs.json').read_text())
    jobs=prepare_jobs(cfg,tables,functions);require(shutil.which('Rscript'),'Rscript required for PCoA/PERMANOVA. Descriptive distance outputs have been saved; install R + vegan + jsonlite and resume.')
    cmd=['Rscript',str(Path(__file__).with_name('statistics.R')),'--jobs',str(root/'jobs.json'),'--out',str(root)]
    with (root/'R_statistics.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    commit(root,signature,[p for p in root.rglob('*') if p.is_file() and p.name!='COMPLETE.json'],{'jobs':len(jobs)});return jobs
