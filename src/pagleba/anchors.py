"""Explicit linked anchor/target distances; no invented independent replicate tests."""
from pathlib import Path
import collections
import numpy as np
from scipy.spatial.distance import braycurtis
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import read_tsv,write_tsv,require,atomic_json,sha
from .statistics import block_matrix

def run(tables,traits,units,out,threshold=.0001,prevalence=.1):
    pairs=read_tsv(units);require(pairs and {'unit','anchor_sample','target_sample','target_type'}<=set(pairs[0]),'Anchor TSV requires unit,anchor_sample,target_sample,target_type');out=Path(out);out.mkdir(parents=True,exist_ok=True);allrows=[]
    for name,path,kind in [(n,p,'taxonomy') for n,p in tables]+[(n,p,'traits') for n,p in traits]:
        ids,features,X,labels=block_matrix(path,kind);keep=np.mean(np.nan_to_num(X)>=threshold,axis=0)>=prevalence if kind=='taxonomy' else np.any(np.nan_to_num(X)>0,axis=0);X=X[:,keep]
        for p in pairs:
            require(p['anchor_sample'] in ids and p['target_sample'] in ids,'Anchor/target sample absent in block '+name);a=X[ids.index(p['anchor_sample'])];b=X[ids.index(p['target_sample'])];bc=float(braycurtis(a,b)) if np.all(np.isfinite(a+b)) and sum(a+b)>0 else 'NA';native=bc if kind=='taxonomy' else float(np.linalg.norm(np.sqrt(a)-np.sqrt(b)));allrows.append(dict(Block=name,**p,Bray_Curtis=bc,Native_Distance=native,Interpretation='Descriptive linked pairs; shared anchors are not independent replicates'))
        rr=[r for r in allrows if r['Block']==name];levels=list(dict.fromkeys(r['target_type'] for r in rr));fig,ax=plt.subplots(figsize=(12,8));fig.subplots_adjust(bottom=.23,right=.75)
        for unit in dict.fromkeys(r['unit'] for r in rr):
            take={r['target_type']:r for r in rr if r['unit']==unit};yy=[float(take[g]['Bray_Curtis']) if g in take and take[g]['Bray_Curtis']!='NA' else np.nan for g in levels];ax.plot(range(len(levels)),yy,'o-',label=unit)
        ax.set_xticks(range(len(levels)),levels,rotation=15);ax.set_ylim(0,1);ax.set_ylabel('Bray–Curtis dissimilarity');ax.set_title(name+' | Explicit anchor-to-target comparisons');ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),frameon=False);ax.grid(alpha=.15);fig.text(.08,.06,'Each line is one supplied anchor unit. Distances sharing an anchor are dependent.\nNo inferential P-value is assigned by this descriptive adapter.',fontsize=10)
        for ext in ['pdf','svg','png']:fig.savefig(out/(name+'_anchor_distances.'+ext),dpi=300)
        plt.close(fig)
    write_tsv(out/'all_anchor_distances.tsv',list(allrows[0]),allrows);atomic_json(out/'provenance.json',{'units_sha256':sha(units),'sources':{n:sha(p) for n,p in tables+traits},'relative_threshold':threshold,'prevalence':prevalence,'inferential_tests':False})
