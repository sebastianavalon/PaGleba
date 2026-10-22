import collections,html,json,math,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,FancyBboxPatch,FancyArrowPatch
from matplotlib.backends.backend_pdf import PdfPages
from .common import RANKS,read_tsv,atomic_json,sha,write_tsv
from .config import resolve
from .statistics import block_matrix,selected_meta

plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def number(x):
    try:return float(x)
    except (ValueError,TypeError):return float('nan')
def save(fig,stem,book=None):
    for ext in ['pdf','svg','png']:fig.savefig(str(stem)+'.'+ext,dpi=300)
    if book:book.savefig(fig)
    plt.close(fig)
def matrix_file(path):
    rows=read_tsv(path);labels=[r['Sample'] for r in rows];return labels,np.array([[number(r[s]) for s in labels] for r in rows])
def table(rows):
    if not rows:return '<p>No records.</p>'
    cols=list(rows[0]);return '<div class="scroll"><table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in cols)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(str(r.get(k,'')))+'</td>' for k in cols)+'</tr>' for r in rows)+'</table></div>'

def trees(cfg,tables,figs):
    metadata=read_tsv(resolve(cfg,cfg['metadata']));tasks=json.loads(resolve(cfg,cfg['tasks']).read_text());assets=[]
    for marker,path in tables:
        rows=read_tsv(path)
        if not rows:continue
        samples=[s for s in rows[0] if s not in ['Feature']+RANKS];X=np.array([[number(r[s]) for s in samples] for r in rows]);totals=X.sum(axis=0);X=np.divide(X,totals[None,:],out=np.zeros_like(X),where=totals[None,:]>0)
        for task in tasks:
            if task.get('blocks') and marker not in task['blocks']:continue
            md=[r for r in selected_meta(task,metadata) if r['sample_id'] in samples];groups=sorted({r[task['group']] for r in md})
            if not groups:continue
            indices={g:[samples.index(r['sample_id']) for r in md if r[task['group']]==g] for g in groups};color={g:plt.get_cmap('tab10')(i%10) for i,g in enumerate(groups)};vectors={};families=set()
            ranks=['Phylum','Class','Order','Family']
            for i,r in enumerate(rows):
                parts=[]
                for rank in ranks:
                    label=r.get(rank,'');parts.append(label if label and not label.startswith('Unclassified') else 'Unresolved '+rank.lower());key=tuple(parts);vectors[key]=vectors.get(key,np.zeros(len(samples)))+X[i]
                families.add(tuple(parts))
            # Filter leaves by the globally specified relative/prevalence threshold; each ancestor vector remains from the full table.
            leaves=sorted(f for f in families if np.mean(vectors[f]>=cfg['statistics']['relative_threshold'])>=cfg['statistics']['prevalence'])
            with PdfPages(figs/f'{task["id"]}__{marker}__family_tree_atlas.pdf') as book:
                for page,start in enumerate(range(0,len(leaves),24),1):
                    leaf=leaves[start:start+24];n=len(leaf);fig,ax=plt.subplots(figsize=(15,15),subplot_kw={'projection':'polar'});fig.subplots_adjust(left=.12,right=.88,top=.88,bottom=.17);ax.set_theta_offset(np.pi/2);ax.set_theta_direction(-1);ax.set_ylim(0,6.8);ax.axis('off');angles=np.linspace(0,2*np.pi,n,endpoint=False);positions={}
                    for depth in range(1,5):
                        for prefix in sorted({f[:depth] for f in leaf}):
                            take=[j for j,f in enumerate(leaf) if f[:depth]==prefix];positions[prefix]=float(np.mean(angles[take]));theta=positions[prefix];r=depth;v=vectors[prefix];means=[np.mean(v[indices[g]]) for g in groups];maxima=[j for j,x in enumerate(means) if np.isclose(x,max(means),rtol=1e-9,atol=1e-12)];c=color[groups[maxima[0]]] if len(maxima)==1 else '#808080';ax.scatter([theta],[r],s=32,c=[c],edgecolors='white',linewidths=.5,zorder=4)
                            if depth>1:
                                parent=positions[prefix[:-1]];arc=np.linspace(parent,theta,40);ax.plot(arc,np.full_like(arc,r-1),color='#a9b5c0',lw=.9);ax.plot([theta,theta],[r-1,r],color='#a9b5c0',lw=.9)
                    phyla=sorted({f[0] for f in leaf});phylum_colors={p:plt.get_cmap('Set3')(j%12) for j,p in enumerate(phyla)}
                    for j,f in enumerate(leaf):
                        theta=angles[j];ax.bar(theta,.17,width=2*np.pi/n*.92,bottom=4.25,color=phylum_colors[f[0]],edgecolor='white');deg=np.degrees(theta);rotation=90-deg;align='left'
                        if 180<deg<360:rotation+=180;align='right'
                        label=f[-1] if not f[-1].startswith('Unresolved') else f[-2]+' / unresolved family';ax.text(theta,4.55,label,ha=align,va='center',rotation=rotation,rotation_mode='anchor',fontsize=8)
                    group_handles=[Patch(color=color[g],label=g+' (n='+str(len(indices[g]))+')') for g in groups]+[Patch(color='#808080',label='Tied maximum')];phy_handles=[Patch(color=phylum_colors[p],label=p) for p in phyla];fig.legend(handles=group_handles,loc='lower left',bbox_to_anchor=(.04,.06),ncol=2,fontsize=9,frameon=False,title='Node colour: greatest group mean relative abundance');fig.legend(handles=phy_handles,loc='lower right',bbox_to_anchor=(.97,.06),ncol=2,fontsize=9,frameon=False,title='Outer band: phylum');fig.suptitle(f'{task["id"]} | {marker} | Taxonomic hierarchy to family | page {page}',fontsize=16);fig.text(.05,.025,'Phylum → Class → Order → Family. This is a taxonomic hierarchy, not a sequence-based phylogeny.\nNode colours are descriptive maxima, not significance tests. Ancestor maxima use all input descendants before leaf filtering.',fontsize=10)
                    stem=figs/f'{task["id"]}__{marker}__family_tree_{page:02d}';save(fig,stem,book);assets.append((task['id'],marker,'Family taxonomic hierarchy',stem.name+'.svg'))
    return assets

def build(cfg,tables,functions,include_trees=False):
    out=Path(cfg['_out']);root=out/'report';figs=root/'figures';figs.mkdir(parents=True,exist_ok=True);assets=[];jobs=json.loads((out/'statistics/jobs.json').read_text()) if (out/'statistics/jobs.json').exists() else []
    for job in jobs:
        if job.get('status_only'):continue
        folder=Path(job['folder']);code=job['task'];block=job['block'];md=read_tsv(job['metadata']);ids,D=matrix_file(folder/'Bray_Curtis.tsv');groups=[r[job['group']] for r in md];levels=sorted(set(groups));palette={g:plt.get_cmap('tab10')(i%10) for i,g in enumerate(levels)};fig,axs=plt.subplots(1,3,figsize=(19,8));fig.subplots_adjust(left=.09,right=.93,bottom=.29,top=.82,wspace=.48)
        for ax,file,title in zip(axs,['Bray_Curtis.tsv','mean_pairwise_Bray_Curtis.tsv','group_mean_profile_Bray_Curtis.tsv'],['Sample-to-sample dissimilarity','Mean pairwise dissimilarity','Group mean-profile dissimilarity']):
            labels,M=matrix_file(folder/file);im=ax.imshow(np.ma.masked_invalid(M),cmap='viridis',vmin=0,vmax=1);ax.set_xticks(range(len(labels)),labels,rotation=90,fontsize=8);ax.set_yticks(range(len(labels)),labels,fontsize=8);ax.set_title(title,fontsize=11)
            if len(labels)<=16:
                for i in range(len(labels)):
                    for j in range(len(labels)):
                        v=M[i,j];ax.text(j,i,f'{v:.2f}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=8,color='white' if np.isfinite(v) and v<.5 else 'black')
        fig.colorbar(im,ax=axs,shrink=.65,pad=.02,label='Bray–Curtis dissimilarity (0–1)');fig.suptitle(f'{code} | {block} | All group-pair comparisons',fontsize=17);fig.text(.05,.05,'Descriptive distances: pairs sharing samples are dependent. Within-group means are undefined for n = 1.\nNative functional-block inference uses square-root Euclidean distances; Bray–Curtis for traits is supplementary.\nGroup sizes: '+', '.join(g+' n='+str(groups.count(g)) for g in levels),fontsize=10);stem=figs/f'{code}__{block}__Bray_Curtis';save(fig,stem);assets.append((code,block,'All pairwise dissimilarities',stem.name+'.svg'))
        if (folder/'alpha_observed_counts.tsv').exists():
            aa=read_tsv(folder/'alpha_observed_counts.tsv');rr=read_tsv(folder/'alpha_rarefied.tsv');fig,axes=plt.subplots(2,3,figsize=(16,10));fig.subplots_adjust(bottom=.16,top=.87,hspace=.5,wspace=.3)
            for rowindex,(data,metrics,titles) in enumerate([(aa,['Observed_Feature_Rows','Shannon','Simpson'],['Observed feature rows','Shannon diversity','Simpson diversity']),(rr,['Observed_Mean','Shannon_Mean','Simpson_Mean'],['Subsampled feature rows','Subsampled Shannon','Subsampled Simpson'])]):
                for ax,metric,title in zip(axes[rowindex],metrics,titles):
                    for gi,g in enumerate(levels):
                        vals=[number(r[metric]) for r in data if r['Group']==g];vals=[v for v in vals if np.isfinite(v)];jitter=np.linspace(-.10,.10,len(vals)) if vals else [];ax.scatter(gi+np.array(jitter),vals,s=45,color=palette[g],edgecolor='white',linewidth=.5)
                    ax.set_xticks(range(len(levels)),[g+'\n(n='+str(groups.count(g))+')' for g in levels],rotation=15,fontsize=9);ax.set_title(title,fontsize=11);ax.grid(axis='y',alpha=.15)
            fig.suptitle(f'{code} | {block} | Alpha diversity per biological sample',fontsize=16);fig.text(.06,.045,'Observed counts are not artificially scaled. Subsampling means come from repeated draws without replacement.\nEach dot is one biological sample; repeated draws are not replicates. Samples below the target depth have no subsampled dot.\nObserved feature rows include unresolved taxonomic bins and are not species or ASVs.',fontsize=10);stem=figs/f'{code}__{block}__alpha';save(fig,stem);assets.append((code,block,'Alpha diversity',stem.name+'.svg'))
        if (folder/'PCoA.tsv').exists():
            xy=read_tsv(folder/'PCoA.tsv');eig=read_tsv(folder/'PCoA_eigenvalues.tsv');fig,ax=plt.subplots(figsize=(12,9));fig.subplots_adjust(right=.73,bottom=.16)
            for g in levels:
                rr=[r for r in xy if r['Group']==g];ax.scatter([number(r['PCoA1']) for r in rr],[number(r['PCoA2']) for r in rr],color=palette[g],s=75,label=g+' (n='+str(len(rr))+')',edgecolor='white')
            for r in xy:ax.annotate(r['Sample'],(number(r['PCoA1']),number(r['PCoA2'])),xytext=(4,5),textcoords='offset points',fontsize=8)
            vals=[number(r['Positive_fraction']) for r in eig if number(r['Eigenvalue'])>0];ax.set_xlabel(f'PCoA 1 ({100*vals[0]:.1f}%)' if vals else 'PCoA 1');ax.set_ylabel(f'PCoA 2 ({100*vals[1]:.1f}%)' if len(vals)>1 else 'PCoA 2');ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),frameon=False);ax.grid(alpha=.15);ax.margins(.2);ax.set_title(f'{code} | {block} | Principal coordinates');fig.text(.08,.04,'Lingoes correction; percentages use positive eigenvalues. Ordination is not an inferential test.',fontsize=10);stem=figs/f'{code}__{block}__PCoA';save(fig,stem);assets.append((code,block,'Principal coordinates',stem.name+'.svg'))
    for taskdir in sorted((out/'models').glob('*')) if (out/'models').exists() else []:
        if not taskdir.is_dir():continue
        code=taskdir.name
        if (taskdir/'network_layout.tsv').exists():
            nd=read_tsv(taskdir/'network_layout.tsv');edges=read_tsv(taskdir/'network_selected_edges.tsv');nodeinfo={r['ID']:r for r in read_tsv(taskdir/'network_nodes.tsv')};xy={r['ID']:(number(r['X']),number(r['Y'])) for r in nd};blocks=sorted({r['Block'] for r in nd});colors={b:plt.get_cmap('tab10')(i%10) for i,b in enumerate(blocks)};fig,ax=plt.subplots(figsize=(14,11));fig.subplots_adjust(bottom=.19)
            for e in edges:
                a,b=xy[e['Source']],xy[e['Target']];ax.plot([a[0],b[0]],[a[1],b[1]],color='#c0392b' if number(e['Rho'])>0 else '#2471a3',alpha=.45,lw=1+2*abs(number(e['Rho'])))
            for r in nd:
                x,y=xy[r['ID']];info=nodeinfo[r['ID']];ax.scatter(x,y,s=100+35*number(r['Degree']),color=colors[r['Block']],marker='o' if info['Kind']=='taxonomy' else 's',edgecolors='white',zorder=3);ax.annotate(r['ID'],(x,y),xytext=(5,5),textcoords='offset points',fontsize=9)
            ax.axis('off');ax.margins(.2);ax.set_title(f'{code} | Integrated association network',fontsize=17);fig.legend(handles=[Patch(color=colors[b],label=b) for b in blocks],loc='lower center',ncol=4,frameon=False);fig.text(.06,.07,'Red: positive; blue: negative partial rank correlation. Node IDs refer to network_nodes.tsv.\nExploratory associations; taxonomically inferred functional traits are not independent measurements. No causal interpretation.',fontsize=10);stem=figs/f'{code}__network';save(fig,stem);assets.append((code,'integrated','Association network',stem.name+'.svg'))
        for pf in sorted(taskdir.glob('PLS_*')):
            if not pf.is_dir() or not (pf/'paths.tsv').exists():continue
            paths=read_tsv(pf/'paths.tsv');n=len(paths);fig,ax=plt.subplots(figsize=(13,max(7,.55*n+3)));fig.subplots_adjust(left=.32,bottom=.19,top=.88)
            for i,r in enumerate(paths):
                x=number(r['Estimate']);y=n-1-i;lo=number(r['CI_low']);hi=number(r['CI_high']);ax.scatter(x,y,s=65,color='#154360',zorder=4)
                if np.isfinite(lo):ax.plot([number(r['Simultaneous_low']),number(r['Simultaneous_high'])],[y,y],color='#aab7b8',lw=1);ax.plot([lo,hi],[y,y],color='#3498db',lw=3)
                ax.annotate(f'{x:.2f}',(x,y),xytext=(0,9),textcoords='offset points',ha='center')
            ax.set_yticks(range(n),[r['Path'] for r in paths][::-1]);ax.axvline(0,color='#85929e',ls='--');ax.set_xlabel('Standardized path coefficient');ax.set_title(f'{code} | PLS-PM ({pf.name[4:]})',loc='left',fontsize=16);ax.margins(y=.2);fig.text(.06,.045,'One local PC per formative block. Blue: 95% conditional bootstrap CI; grey: simultaneous max-|t| CI.\nNo interval is drawn for descriptive-only estimates. Paths are model assumptions; taxonomy-derived functions do not prove causation.',fontsize=10);stem=figs/f'{code}__{pf.name}';save(fig,stem);assets.append((code,'integrated','PLS-PM coefficients',stem.name+'.svg'))
            blocks=list(dict.fromkeys(v for r in paths for v in r['Path'].split(' -> ')));depth={b:0 for b in blocks};
            for _ in blocks:
                for r in paths:
                    a,b=r['Path'].split(' -> ');depth[b]=max(depth[b],depth[a]+1)
            layers={d:[b for b in blocks if depth[b]==d] for d in set(depth.values())};height=max(1,max(len(v) for v in layers.values())-1);pos={b:(3*d,float(y)) for d,bs in layers.items() for b,y in zip(bs,np.linspace(0,height,len(bs)))};fig,ax=plt.subplots(figsize=(max(15,4+5*max(depth.values())),9));ax.set(xlim=(-.8,3*max(depth.values())+.8),ylim=(-.8,height+.8));ax.axis('off')
            for i,r in enumerate(paths):
                a,b=r['Path'].split(' -> ');x,y=pos[a];xx,yy=pos[b];sig=number(r['Simultaneous_low'])*number(r['Simultaneous_high'])>0;ax.add_patch(FancyArrowPatch((x+.55,y),(xx-.55,yy),arrowstyle='-|>',mutation_scale=15,color='#1a5276' if sig else '#85929e',lw=2,linestyle='-' if sig else '--'));fraction=.35 if i%2 else .65;px=(x+.55)+(xx-x-1.10)*fraction;py=y+(yy-y)*fraction;ax.text(px,py+.09,f'β={number(r["Estimate"]):.2f}',ha='center',bbox=dict(facecolor='white',edgecolor='none'),fontsize=11)
            for i,(b,(x,y)) in enumerate(pos.items()):ax.add_patch(FancyBboxPatch((x-.5,y-.15),1,.3,boxstyle='round,pad=.04',facecolor=plt.get_cmap('tab10')(i%10),edgecolor='none'));ax.text(x,y,b,ha='center',va='center',color='white',fontsize=11)
            ax.set_title(f'{code} | PLS-PM paths ({pf.name[4:]})',fontsize=17);fig.text(.06,.04,'Solid arrows: simultaneous CI excludes zero. Dashed: not confirmed or descriptive only.\nArrows describe a prespecified statistical model, not demonstrated causation.',fontsize=11);stem=figs/f'{code}__{pf.name}__paths';save(fig,stem);assets.append((code,'integrated','PLS-PM path diagram',stem.name+'.svg'))
    for hitfile in sorted((out/'LEfSe').rglob('biomarkers.tsv')) if (out/'LEfSe').exists() else []:
        hits=read_tsv(hitfile);hits=sorted(hits,key=lambda r:-abs(number(r['LDA_Log10'])));task=hitfile.parent.parent.name;block=hitfile.parent.name;groups=sorted({r['Enriched_Group'] for r in hits});colors={g:plt.get_cmap('tab10')(i%10) for i,g in enumerate(groups)}
        for start in range(0,len(hits),28):
            take=hits[start:start+28];fig,ax=plt.subplots(figsize=(15,max(7,.35*len(take)+3)));fig.subplots_adjust(left=.43,bottom=.15,top=.88);y=np.arange(len(take));ax.barh(y,[abs(number(r['LDA_Log10'])) for r in take],color=[colors[r['Enriched_Group']] for r in take]);ax.set_yticks(y,[r['Taxonomy'].split('|')[-1][3:]+(' *' if r['BH_Supported']=='True' else '') for r in take],fontsize=9);ax.invert_yaxis();ax.set_xlabel('Official LEfSe LDA score (log10)');ax.set_title(f'{task} | {block} | LEfSe biomarkers',fontsize=16);fig.legend(handles=[Patch(color=colors[g],label=g) for g in groups],loc='lower center',ncol=4,frameon=False);fig.text(.06,.06,'* Kruskal–Wallis BH q ≤ 0.05 across all tested features. Classical LEfSe hits and FDR support are reported separately.\nLEfSe subclass consistency is not equivalent to adjustment by a covariate regression model.',fontsize=10);stem=figs/f'{task}__{block}__LEfSe_{start//28+1:02d}';save(fig,stem);assets.append((task,block,'Official LEfSe biomarkers',stem.name+'.svg'))
    if include_trees:assets.extend(trees(cfg,tables,figs))
    stat=read_tsv(out/'statistics/PERMANOVA_summary.tsv') if (out/'statistics/PERMANOVA_summary.tsv').exists() else [];model=read_tsv(out/'models/MODEL_STATUS.tsv') if (out/'models/MODEL_STATUS.tsv').exists() else []
    qcs=[(m,json.loads((p.parent/'count_reconciliation.json').read_text())) for m,p in tables if (p.parent/'count_reconciliation.json').exists()]
    sections=''.join('<section><h2>'+html.escape(t+' | '+b+' | '+title)+'</h2><object type="image/svg+xml" data="figures/'+f+'"></object><p><a href="figures/'+f[:-4]+'.pdf">Vector PDF</a> · <a href="figures/'+f[:-4]+'.png">PNG 300 dpi</a></p></section>' for t,b,title,f in assets)
    links=''.join('<li><a href="../'+html.escape(str(p.relative_to(out)))+'">'+html.escape(str(p.relative_to(out)))+'</a></li>' for p in sorted(out.rglob('*.tsv')) if 'report' not in p.relative_to(out).parts)
    text='<!doctype html><html lang="en"><meta charset="utf-8"><title>PaGleba analysis report</title><style>body{font:16px Arial;background:#f1f5f8;color:#17324d;max-width:1450px;margin:auto;padding:30px}section{background:white;padding:20px;margin:25px 0;border-radius:12px}object{width:100%;height:800px}.scroll{overflow:auto}table{border-collapse:collapse}td,th{padding:7px;border:1px solid #ccd5dd;font-size:13px}a{color:#176ba0}</style><h1>PaGleba reproducible ONT analysis</h1><p>Database release names are retained. The pipeline does not claim that all nomenclature has independently been updated to current registries. Taxonomic tables are read-occurrence classifications, not ASV or sequence-cluster OTUs. All source reads and explicit technical/assignment statuses are accounted for.</p><p>Taxonomic exclusions: '+html.escape(json.dumps(cfg.get('taxonomy_filter',{}).get('exclude_clades',[])))+'. Plant reads are retained unless an explicit filter is configured. Functional predictions describe reference-derived potential, not measured activity. Small/unreplicated comparisons are descriptive; pairwise distances are not independent replicates.</p><h2>PERMANOVA status and results</h2>'+table(stat)+'<h2>Task-specific model status</h2>'+table(model)+'<h2>Count reconciliation</h2><pre>'+html.escape(json.dumps(dict(qcs),indent=2))+'</pre><h2>Tables</h2><ul>'+links+'</ul>'+sections+'</html>'
    if cfg.get('input_tables'):text=text.replace('<h2>PERMANOVA status and results</h2>','<p>Imported tables are analyzed as supplied. Upstream host/taxon filtering and artificial rescaling cannot be inferred from the matrix; their provenance must be documented separately.</p><h2>PERMANOVA status and results</h2>')
    if cfg.get('demonstration_only'):text=text.replace('<h1>PaGleba reproducible ONT analysis</h1>','<h1>PaGleba synthetic software demonstration</h1><p><strong>SYNTHETIC TEST DATA: this is not a biological study or an assay validation.</strong></p>')
    (root/'REPORT.html').write_text(text);atomic_json(root/'assets.json',assets);atomic_json(root/'validation.json',{'missing_assets':[f for t,b,title,f in assets if not (figs/f).exists()],'figures':len(assets),'language':'English','report_sha256':sha(root/'REPORT.html')});return root/'REPORT.html'
