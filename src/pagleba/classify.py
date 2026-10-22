"""ONT minimap2 + conservative competing-hit LCA, one auditable record per read occurrence."""
import collections,functools,gzip,itertools,json,math,os,signal,sqlite3,subprocess,time
from pathlib import Path
from .common import RANKS,WorkflowError,require,sha,digest,atomic_json,checked_checkpoint,commit,fastq,write_tsv
from .config import resolve

FIELDS=['QueryID','ReadID','Sample','Length','Status','AssignedRank','LevelsJSON','Identity','QueryCoverage','AlignmentScore','BestReference','CompetitiveHits','CandidateLimitReached']
def acceptable_genus(s):
    return bool(s) and not any(v in s.casefold() for v in ['unclassified','incertae','uncultured','unidentified','unknown','.gen'])
def parse_paf(line):
    f=line.rstrip('\n').split('\t');require(len(f)>=12,'Invalid PAF');tags={v[:2]:v[5:] for v in f[12:]}
    return {'query':f[0],'ref':f[5],'score':int(tags.get('AS',0)),'identity':int(f[9])/int(f[10]) if int(f[10]) else 0,'coverage':(int(f[3])-int(f[2]))/int(f[1]) if int(f[1]) else 0}
def common_levels(paths):
    result=[]
    for values in zip(*paths):
        if len(set(values))!=1:break
        result.append(values[0])
    return result[:7]+['']*max(0,8-len(result[:7]))
def classify_hits(hits,get_levels,m):
    valid=[h for h in hits if h['identity']>=m['identity'] and h['coverage']>=m['coverage']]
    best=max(valid or hits,key=lambda h:(h['score'],h['identity'],h['coverage']))
    if not valid:return {'Status':'no_acceptable_alignment','LevelsJSON':'[]','AssignedRank':'','Identity':best['identity'],'QueryCoverage':best['coverage'],'AlignmentScore':best['score'],'BestReference':best['ref'],'CompetitiveHits':0,'CandidateLimitReached':False}
    cutoff=best['score']*m['score_ratio'];top=[h for h in valid if h['score']>=cutoff];levels=common_levels([get_levels(h['ref']) for h in top]);limit=len(hits)>=m.get('secondary',50)+1 and min(h['score'] for h in hits)>=cutoff
    if limit:levels[6:]=['',''] # Withhold genus if the near-best candidate list is truncated.
    names=';'.join(levels).casefold();domain=levels[0]
    if any(v in names for v in ['_mitochondrion','_plastid','_apicoplast','_nucleomorph','_chromatophore','chloroplast','mitochondria']):status='organelle'
    elif '_pseudogene' in names:status='pseudogene'
    elif 'unispike' in names:status='synthetic_control'
    elif domain=='Eukaryota':status='eukaryotic_genus' if acceptable_genus(levels[6]) else 'eukaryotic_unresolved_genus'
    elif domain in {'Bacteria','Archaea'}:status=domain.casefold()+'_genus' if acceptable_genus(levels[6]) else domain.casefold()+'_unresolved_genus'
    else:status='unresolved_domain'
    rank=next((RANKS[j] for j in range(6,-1,-1) if levels[j]),'')
    return {'Status':status,'LevelsJSON':json.dumps(levels),'AssignedRank':rank,'Identity':best['identity'],'QueryCoverage':best['coverage'],'AlignmentScore':best['score'],'BestReference':best['ref'],'CompetitiveHits':len(top),'CandidateLimitReached':limit}

def mean_qscore(q):return -10*math.log10(sum(10**(-(ord(c)-33)/10) for c in q)/len(q))

def classify_chunk(folder,records,sample,start,m,db,dbmeta,threads):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);signature=digest({'algorithm':'minimap2-lca-v1','code_sha256':sha(Path(__file__)),'database':dbmeta,'records':[(h,s,q) for h,s,p,q in records],'sample':sample,'start':start,'parameters':m})
    if checked_checkpoint(folder,signature):return {'resumed':True,'reads':len(records)}
    sql=sqlite3.connect('file:'+str(db/'taxonomy.sqlite')+'?mode=ro',uri=True)
    @functools.lru_cache(maxsize=100000)
    def get_levels(ref):
        row=sql.execute('SELECT levels FROM refs WHERE reference=?',(ref,)).fetchone();require(row is not None,f'Alignment reference absent in taxonomy: {ref}');return json.loads(row[0])
    meta={};initial=[];query=folder/'queries.fasta';query_count=0
    with query.open('w') as f:
        for i,(h,s,p,q) in enumerate(records,start):
            qid=f'{sample}_r{i:08d}';meta[qid]={'QueryID':qid,'ReadID':h[1:].split()[0],'Sample':sample,'Length':len(s)}
            status='length_filtered' if not m['min_length']<=len(s)<=m['max_length'] else 'quality_filtered' if m.get('min_read_qscore') is not None and mean_qscore(q)<m['min_read_qscore'] else None
            if status:initial.append(dict(meta[qid],Status=status,LevelsJSON='[]'))
            else:f.write('>'+qid+'\n'+s+'\n');query_count+=1
    cmd=['minimap2','-x',m.get('preset','map-ont'),'-t',str(threads),'-K','100M','-c','--secondary=yes','-N',str(m.get('secondary',50)),'-p','0.90',str(db/'reference.mmi'),str(query)];atomic_json(folder/'command.json',cmd)
    reported=set();statuses=collections.Counter();assign=folder/'assignments.tsv.gz';tmp=folder/'assignments.partial.tsv.gz'
    import csv
    with gzip.open(tmp,'wt',compresslevel=1,newline='') as fa,gzip.open(folder/'alignments.paf.gz','wt',compresslevel=1) as paf,(folder/'minimap2.log').open('w') as log:
        w=csv.DictWriter(fa,FIELDS,delimiter='\t',extrasaction='ignore');w.writeheader()
        def emit(row):
            qid=row['QueryID'];require(qid in meta and qid not in reported,f'Unexpected/repeated query ID {qid}; multi-part indexes are not supported here');reported.add(qid);row.setdefault('CandidateLimitReached',False);row.setdefault('CompetitiveHits',0);w.writerow(row);statuses[row['Status']]+=1
        for row in initial:emit(row)
        if query_count:
            process=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log,text=True,start_new_session=True)
            try:
                def stream():
                    for line in process.stdout:paf.write(line);yield parse_paf(line)
                for qid,hitlist in itertools.groupby(stream(),key=lambda h:h['query']):
                    hits=list(hitlist);require(qid in meta,f'Unknown mapper query {qid}');emit(dict(meta[qid],**classify_hits(hits,get_levels,m)))
                require(process.wait()==0,f'minimap2 failed: {folder/"minimap2.log"}')
            except BaseException:
                if process.poll() is None:
                    os.killpg(process.pid,signal.SIGTERM)
                    try:process.wait(timeout=10)
                    except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                raise
            finally:
                process.stdout.close()
        for qid,row in meta.items():
            if qid not in reported:emit(dict(row,Status='unaligned',LevelsJSON='[]'))
    sql.close();require(len(reported)==len(records),'Read count mismatch');tmp.replace(assign);query.unlink(missing_ok=True)
    # Read gzip files to EOF before marking atomic completion (CRC validation).
    for f in [assign,folder/'alignments.paf.gz']:
        with gzip.open(f,'rb') as z:
            for _ in iter(lambda:z.read(8*1024*1024),b''):pass
    atomic_json(folder/'qc.json',{'reads':len(records),'statuses':dict(statuses),'eligible_for_alignment':query_count});commit(folder,signature,[assign,folder/'alignments.paf.gz',folder/'qc.json',folder/'command.json',folder/'minimap2.log'])
    return {'resumed':False,'reads':len(records)}

def classify(cfg,trimmed,pilot=None):
    results=[];start_time=time.monotonic();new_reads=0;new_chunks=0;out=Path(cfg['_out'])/'taxonomy';versions={};verified_db={}
    for item in trimmed:
        marker=item['marker'];m=cfg['markers'][marker];db=resolve(cfg,m['database'])
        if str(db) not in verified_db:
            meta=json.loads((db/'database.json').read_text());require(meta['marker']==m['marker'],f'Marker mismatch: {marker} cannot use {meta["marker"]} database')
            require(meta['preset']==m.get('preset','map-ont'),'Index preset differs from mapper preset; rebuild the database index')
            require(meta['provider']==m['provider'],'Database provider mismatch')
            require(meta['version']==m['version'],'Requested database release differs from prepared database metadata')
            for f,h in meta['files'].items():require(sha(db/f)==h,f'Database file changed: {db/f}')
            actual=subprocess.check_output(['minimap2','--version'],text=True).strip();require(actual==meta['minimap2_version'],'minimap2 version changed since indexing; rebuild in a new database folder')
            verified_db[str(db)]=meta
        meta=verified_db[str(db)];source=Path(item['fastq']);sample=item['sample_id'];folder=out/marker/sample;folder.mkdir(parents=True,exist_ok=True)
        config={'input_sha256':sha(source),'parameters':m,'chunk_reads':cfg['chunk_reads'],'database_sha256':sha(db/'database.json')}
        prev=folder/'configuration.json'
        if prev.exists():require(json.loads(prev.read_text())==config,f'Changed sample taxonomy configuration: {sample}; use a new output project')
        else:atomic_json(prev,config)
        it=iter(fastq(source));start=1;chunks=[]
        while True:
            batch=list(itertools.islice(it,cfg['chunk_reads']))
            if not batch:break
            ck=folder/'checkpoints'/f'{start:012d}_{start+len(batch)-1:012d}';r=classify_chunk(ck,batch,sample,start,m,db,meta,cfg['threads']);chunks.append(str(ck));start+=len(batch)
            if not r['resumed']:new_chunks+=1;new_reads+=r['reads']
            progress={'sample':sample,'marker':marker,'completed_new_chunks':new_chunks,'new_reads':new_reads,'elapsed_seconds':round(time.monotonic()-start_time,2),'total_trimmed_reads':sum(v['reads'] for v in trimmed),'note':'Pilot duration includes input/database verification and is only an estimate; samples differ.'};atomic_json(out/'progress.json',progress);print('CHECKPOINT',marker,sample,ck.name,flush=True)
            if pilot and new_chunks>=pilot:
                atomic_json(out/'pilot.json',progress);return None
        require(start-1==item['reads'],'FASTQ count changed during classification');results.append(dict(item,chunks=chunks))
    atomic_json(out/'manifest.json',results);atomic_json(out/'CLASSIFICATION_COMPLETE.json',{'samples':len(results),'read_occurrences':sum(v['reads'] for v in results)});return results
