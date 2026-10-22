import collections,csv,gzip,json
from pathlib import Path
from .common import RANKS,require,sha,digest,atomic_json,write_tsv,commit,checked_checkpoint

def export_tables(cfg,classified):
    outputs=[]
    for marker in sorted({v['marker'] for v in classified}):
        m=cfg['markers'][marker];items=[v for v in classified if v['marker']==marker];names=[v['sample_id'] for v in items];require(len(names)==len(set(names)),'Repeated marker/sample output');folder=Path(cfg['_out'])/'tables'/marker;folder.mkdir(parents=True,exist_ok=True)
        chunks=[Path(c) for v in items for c in v['chunks']];sig=digest({'algorithm':'aggregate-v1','code_sha256':sha(Path(__file__)),'chunks':{str(p):sha(p/'COMPLETE.json') for p in chunks},'samples':names,'filter':cfg.get('taxonomy_filter',{}),'domains':m.get('domains')})
        if checked_checkpoint(folder,sig):outputs.append((marker,folder/'analysis_counts.tsv'));continue
        bins=collections.defaultdict(collections.Counter);status_counts=collections.defaultdict(collections.Counter);seen=collections.Counter();limits=collections.Counter();allpairs=[]
        for item in items:
            sample=item['sample_id'];expected=0
            for ck in item['chunks']:
                ck=Path(ck)
                with gzip.open(ck/'assignments.tsv.gz','rt') as f:
                    for r in csv.DictReader(f,delimiter='\t'):
                        seen[sample]+=1;expected+=1;status_counts[r['Status']][sample]+=1;levels=json.loads(r['LevelsJSON'] or '[]');limits[sample]+=r.get('CandidateLimitReached')=='True'
                        if levels:bins[(r['Status'],tuple(levels))][sample]+=1
            require(expected==item['reads'],f'Count reconciliation failed for {sample}')
        domains=m.get('domains') or (['Bacteria'] if m['marker']=='16s' else ['Eukaryota'])
        rows=[];mapping=[];alltax=[]
        for (status,levels),counts in sorted(bins.items()):
            feature='PG_'+digest({'status':status,'levels':levels})[:16];fields=dict(zip(RANKS,(list(levels)+['']*8)[:8]));fields['Species']=''
            if status.endswith('unresolved_genus'):fields['Genus']='Unclassified_genus'
            row=dict(Feature=feature,**fields,**{s:counts[s] for s in names});alltax.append(dict(row,Status=status))
            if levels[0] in domains and (status.endswith('_genus')):
                rows.append(row);mapping.append(dict(Feature=feature,Status=status,OriginalLevelsJSON=json.dumps(levels),PathID_Provenance='Local SHA256 lineage hash; not NCBI taxid'))
        fields=['Feature']+RANKS+names;write_tsv(folder/'otu_genus_counts.tsv',fields,rows);write_tsv(folder/'feature_taxonomy.tsv',['Feature','Status','OriginalLevelsJSON','PathID_Provenance'],mapping)
        write_tsv(folder/'taxonomy_all_domains_counts.tsv',['Feature']+RANKS+['Status']+names,alltax)
        write_tsv(folder/'all_status_counts.tsv',['Status']+names,[dict(Status=status,**{s:counts[s] for s in names}) for status,counts in sorted(status_counts.items())])
        excluded={s.casefold() for s in cfg.get('taxonomy_filter',{}).get('exclude_clades',[])};original={r['Feature']:json.loads(r['OriginalLevelsJSON']) for r in mapping}
        dropped=[r for r in rows if any(v.casefold() in excluded for v in original[r['Feature']])];kept=[r for r in rows if r not in dropped]
        write_tsv(folder/'analysis_counts.tsv',fields,kept);write_tsv(folder/'excluded_taxa.tsv',fields,dropped)
        for filename,rs in [('otu_genus_relative_abundance.tsv',rows),('analysis_relative_abundance.tsv',kept)]:
            total={s:sum(r[s] for r in rs) for s in names};relative=[dict(r,**{s:r[s]/total[s] if total[s] else 0 for s in names}) for r in rs];write_tsv(folder/filename,fields,relative)
        statusfields=sorted(status_counts);qc=[]
        for it in items:
            s=it['sample_id'];ret=sum(r[s] for r in rows);an=sum(r[s] for r in kept);require(sum(c[s] for c in status_counts.values())==seen[s]==it['reads'],'Status reconciliation failed')
            qc.append(dict(Sample=s,Merged_Reads=it['input_reads'],After_Primer_QC=it['reads'],Target_Domain_Reads=ret,Analysis_Reads=an,Taxon_Filter_Removed=ret-an,Candidate_Limit_Reached=limits[s],**{k:status_counts[k][s] for k in statusfields}))
        write_tsv(folder/'sample_QC.tsv',list(qc[0]) if qc else ['Sample'],qc);atomic_json(folder/'count_reconciliation.json',{'passed':True,'samples':len(names),'target_domains':domains,'source_reads':sum(seen.values()),'target_table_reads':sum(r[s] for r in rows for s in names),'analysis_reads':sum(r[s] for r in kept for s in names),'excluded_clades':sorted(excluded),'duplicates':'all occurrences retained','singleton_filter':False,'artificial_count_scaling':False,'table_type':'taxonomic abundance; not ASV or sequence-cluster OTUs','empty_target_samples':[q['Sample'] for q in qc if q['Target_Domain_Reads']==0]})
        commit(folder,sig,[f for f in folder.iterdir() if f.is_file() and f.name!='COMPLETE.json']);outputs.append((marker,folder/'analysis_counts.tsv'))
    return outputs
