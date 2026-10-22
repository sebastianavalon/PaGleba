import collections,csv,json,re,shutil,subprocess
from pathlib import Path
import numpy as np
from scipy.stats import kruskal
from .common import RANKS,read_tsv,write_tsv,require,atomic_json,sha,digest,checked_checkpoint,commit
from .statistics import selected_meta
from .config import resolve

def bh(p):
    p=np.array(p,float);order=np.argsort(p);q=np.ones(len(p));last=1.
    for j in range(len(p)-1,-1,-1):i=order[j];last=min(last,p[i]*len(p)/(j+1));q[i]=last
    return q

def run_lefse(cfg,tables):
    spec=cfg.get('lefse',{});root=Path(cfg['_out'])/'LEfSe';root.mkdir(parents=True,exist_ok=True);minimum=spec.get('min_class',5);summaries=[]
    if not spec.get('enabled',False):atomic_json(root/'status.json',{'status':'not_configured','reason':'Enable official LEfSe explicitly; no log-mean surrogate is labelled LDA'});return
    fmt=spec.get('format_executable') or shutil.which('lefse_format_input.py');run=spec.get('run_executable') or shutil.which('lefse_run.py');require(fmt and run,'Official LEfSe executables required; see environment-lefse.yml')
    metadata=read_tsv(resolve(cfg,cfg['metadata']));tasks=json.loads(resolve(cfg,cfg['tasks']).read_text())
    for block,source in tables:
        rows=read_tsv(source);require(rows,'Empty LEfSe input');samples=[k for k in rows[0] if k not in ['Feature']+RANKS];X=np.array([[float(r[s]) for r in rows] for s in samples]);totals=X.sum(axis=1);X=np.divide(X,totals[:,None],out=np.zeros_like(X),where=totals[:,None]>0);agg={}
        for k,r in enumerate(rows):
            path=[]
            for prefix,rank in zip('dkpcofg',RANKS[:-1]):
                v=r.get(rank,'');known=bool(v) and not any(t in v.casefold() for t in ['unclassified','incertae','uncultured','unknown']);path.append(prefix+'__'+(re.sub(r'[^A-Za-z0-9_.-]','_',v) if known else 'unresolved_'+prefix))
                if known:
                    key='|'.join(path);agg[key]=agg.get(key,np.zeros(len(samples)))+X[:,k]
        # Ancestors are summed before any descendant prevalence filter. Constant domains are not tested.
        agg={k:v for k,v in agg.items() if np.mean(v>=cfg['statistics']['relative_threshold'])>=cfg['statistics']['prevalence'] and np.std(v)>1e-12}
        for t in tasks:
            if spec.get('tasks') and t['id'] not in spec['tasks']:continue
            md=selected_meta(t,metadata);folder=root/t['id']/block;folder.mkdir(parents=True,exist_ok=True);g=[r[t['group']] for r in md];counts=collections.Counter(g);status={'Task':t['id'],'Block':block,'N':len(md),'Status':'not_inferred_small_or_missing_group','Classical_hits':0,'BH_supported_hits':0}
            if len(counts)<2 or min(counts.values())<minimum or any(r['sample_id'] not in samples or totals[samples.index(r['sample_id'])]<=0 for r in md):summaries.append(status);continue
            ix=[samples.index(r['sample_id']) for r in md];features=[(k,v[ix]) for k,v in sorted(agg.items()) if np.std(v[ix])>1e-12];require(features,'No nonconstant LEfSe features');sig=digest({'algorithm':'official-lefse-v1','table':sha(source),'metadata':md,'features':[k for k,v in features],'settings':spec})
            if checked_checkpoint(folder,sig):summaries.append(json.loads((folder/'summary.json').read_text()));continue
            with (folder/'input.tsv').open('w') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(['Class']+g);w.writerow(['Subclass']+[r[t['strata']] if t.get('strata') else r[t['group']] for r in md]);w.writerow(['SampleID']+[r['sample_id'] for r in md])
                for i,(tax,v) in enumerate(features):w.writerow([f'T{i:06d}']+list(1e6*v))
            cmds=[[fmt,str(folder/'input.tsv'),str(folder/'input.in'),'-c','1','-s','2','-u','3','-o','-1'],[run,str(folder/'input.in'),str(folder/'lefse.res'),'-a','.05','-w','.05','-l',str(spec.get('lda',2)),'-b','30','-e','1','--min_c',str(minimum),'-y','0']];atomic_json(folder/'commands.json',cmds)
            with (folder/'run.log').open('w') as log:
                for cmd in cmds:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
            p=[]
            for tax,v in features:
                try:p.append(kruskal(*[v[np.array(g)==group] for group in sorted(counts)]).pvalue)
                except ValueError:p.append(1.)
            q=bh(p);out=[]
            for line in (folder/'lefse.res').read_text().splitlines():
                vals=line.split('\t')+['']*5;fid,logmean,group,lda,pval=vals[:5];j=int(fid[1:]);out.append(dict(Feature=fid,Taxonomy=features[j][0],Enriched_Group=group,LDA_Log10=lda,KW_P=p[j],KW_BH_Q=q[j],Classical_Hit=bool(group and lda),BH_Supported=bool(group and lda and q[j]<=.05)))
            require(len(out)==len(features),'Official LEfSe did not return all tested features');write_tsv(folder/'all_features.tsv',list(out[0]),out);write_tsv(folder/'biomarkers.tsv',list(out[0]),[r for r in out if r['Classical_Hit']]);status.update(Status='official_LEfSe_exploratory',Classical_hits=sum(r['Classical_Hit'] for r in out),BH_supported_hits=sum(r['BH_Supported'] for r in out));atomic_json(folder/'summary.json',status);commit(folder,sig,[p for p in folder.iterdir() if p.is_file() and p.name!='COMPLETE.json']);summaries.append(status)
    write_tsv(root/'summary.tsv',['Task','Block','N','Status','Classical_hits','BH_supported_hits'],summaries)
