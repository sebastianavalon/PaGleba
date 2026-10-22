"""Small synthetic project; never presents test data as biological validation."""
import csv,gzip,json,random
from pathlib import Path
from .common import atomic_json,write_tsv,require,RANKS
from .database import prepare
from .config import initialize

def create(folder,n=12,threads=2):
    folder=Path(folder).resolve();require(n>=4 and n%2==0,'Choose an even number of synthetic samples, at least four');cfgpath=initialize(folder,'18s');cfg=json.loads(cfgpath.read_text());cfg['project']='synthetic_software_demo';cfg['threads']=threads;cfg['chunk_reads']=15;cfg['statistics'].update(permutations=99,bootstrap=40,rarefaction_depth=20,rarefaction_repeats=20);cfg['demonstration_only']=True;rows=[];metadata=[];rng=random.Random(87450)
    markers={}
    for marker,length in [('ssu',850),('its',650)]:
        source=folder/'reference_inputs'/marker;source.mkdir(parents=True,exist_ok=True);references={};tax=[]
        for i in range(2):
            name='ref'+str(i+1);references[name]=''.join(rng.choice('ACGT') for _ in range(length));lineage=['Eukaryota','Viridiplantae' if i==0 else 'Fungi','Streptophyta' if i==0 else 'Ascomycota','DemoClass'+str(i),'DemoOrder'+str(i),'DemoFamily'+str(i),'DemoGenus'+str(i),'DemoSpecies'+str(i)];tax.append(dict(reference=name,**dict(zip(RANKS,lineage))))
        (source/'reference.fasta').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in references.items()));write_tsv(source/'taxonomy.tsv',['reference']+RANKS,tax);db=folder/'databases'/marker;prepare(source/'reference.fasta',source/'taxonomy.tsv',db,'eukaryome',marker,'synthetic-demo','tsv',threads=threads)
        m=dict(cfg['markers']['ssu']);m.update(marker=marker,database='databases/'+marker,version='synthetic-demo',min_length=150,max_length=1500);markers[marker]=m
        for j in range(n):
            sample=f'sample{j+1:02d}';group='control' if j<n//2 else 'treatment';p=.25 if group=='control' else .75
            for run in ['run1','run2']:
                path=folder/'raw'/marker/run/sample/'reads.fastq.gz';path.parent.mkdir(parents=True,exist_ok=True)
                with gzip.open(path,'wt') as f:
                    for read in range(12):
                        rid='ref0' if rng.random()<p else 'ref1';seq=references['ref1' if rid=='ref0' else 'ref2'];f.write('@repeated_uuid barcode=barcode01 barcode_alias=barcode01\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
                rows.append(dict(sample_id=sample,marker=marker,run_id=run,barcode='barcode01',fastq=str(path.relative_to(folder)),header_barcode=f'barcode{j+1:02d}'))
    for j in range(n):metadata.append(dict(sample_id=f'sample{j+1:02d}',group='control' if j<n//2 else 'treatment',location='site1' if j%2==0 else 'site2',material='soil',stage='sampling'))
    cfg['markers']=markers;cfg['functions']=[]
    for marker in markers:
        ref=folder/(marker+'_traits.tsv');write_tsv(ref,['Source','Rank','Taxon','Trait','Weight','Evidence','ContextFamily'],[dict(Source='Synthetic',Rank='Genus',Taxon='DemoGenus'+str(i),Trait='example_trait_'+str(i),Weight=1,Evidence='Synthetic software fixture, not biological trait evidence',ContextFamily='') for i in range(2)]);cfg['functions'].append({'id':marker+'_functions','marker':marker,'annotations':ref.name,'method':'long'})
    cfg['models']={'blocks':['ssu','its','ssu_functions','its_functions'],'block_order':['ssu','its','ssu_functions','its_functions'],'paths':[{'from':a,'to':b} for a in ['ssu','its'] for b in ['ssu_functions','its_functions']],'bootstrap_strata':'group'}
    atomic_json(cfgpath,cfg);write_tsv(folder/'samples.tsv',['sample_id','marker','run_id','barcode','fastq','header_barcode'],rows);write_tsv(folder/'metadata.tsv',list(metadata[0]),metadata);atomic_json(folder/'tasks.json',[dict(id='all_groups',group='group',subset={},covariates=[],strata=None,blocks=['ssu','its','ssu_functions','its_functions'],note='Synthetic test only'),dict(id='small_subset',group='group',subset={'sample_id':['sample01',f'sample{n:02d}']},covariates=[],strata=None,blocks=['ssu','its'],note='One sample per group: descriptive only')]);return cfgpath
