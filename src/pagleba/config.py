import json,glob,os,shutil,re
from pathlib import Path
from .common import require,safe_id,read_tsv,atomic_json,write_tsv,WorkflowError

PROFILES={
 '18s':{'marker':'ssu','provider':'eukaryome','version':'2.1','min_length':None,'max_length':None,'identity':.95,'coverage':.8,'score_ratio':.99},
 'its':{'marker':'its','provider':'eukaryome','version':'2.1','min_length':None,'max_length':None,'identity':.95,'coverage':.8,'score_ratio':.99},
 '16s':{'marker':'16s','provider':'silva','version':'144','min_length':1000,'max_length':2000,'identity':.95,'coverage':.8,'score_ratio':.99},
 'lsu':{'marker':'lsu','provider':'eukaryome','version':'2.1','min_length':None,'max_length':None,'identity':.95,'coverage':.8,'score_ratio':.99},
 'long':{'marker':'long','provider':'eukaryome','version':'2.1','min_length':None,'max_length':None,'identity':.95,'coverage':.8,'score_ratio':.99}}

PROFILES['18s-f566-r1200']=dict(PROFILES['18s'],min_length=400,max_length=1200)

def initialize(folder,profile):
    folder=Path(folder).resolve();folder.mkdir(parents=True,exist_ok=True);dest=folder/'project.json';require(not dest.exists(),f'Already exists: {dest}')
    m=dict(PROFILES[profile]);m.update(database='databases/'+m['marker'],trim='none',forward=None,reverse=None,error_rate=.15,overlap=16,preset='map-ont',secondary=50)
    cfg={'schema_version':1,'project':'new_project','samples':'samples.tsv','metadata':'metadata.tsv','tasks':'tasks.json','output':'results','threads':12,'chunk_reads':5000,'seed':20261021,'markers':{m['marker']:m},'taxonomy_filter':{'exclude_clades':[]},'functions':[],'statistics':{'permutations':9999,'relative_threshold':.0001,'prevalence':.10,'network_rho':.6,'network_q':.05,'bootstrap':1000}}
    atomic_json(dest,cfg)
    write_tsv(folder/'samples.tsv',['sample_id','marker','run_id','barcode','fastq','header_barcode'],[{'sample_id':'sample01','marker':m['marker'],'run_id':'run1','barcode':'barcode01','fastq':'/replace/with/run1/barcode01/*.fastq.gz','header_barcode':'barcode01'}])
    write_tsv(folder/'metadata.tsv',['sample_id','group','location','material','stage'],[dict(sample_id='sample01',group='control',location='site1',material='soil',stage='flowering')])
    atomic_json(folder/'tasks.json',[{'id':'all_groups','group':'group','subset':{},'covariates':[],'strata':None,'blocks':[m['marker']],'note':'Replace with the actual experimental question; do not invent replicates.'}])
    return dest

def load(path,threads=None):
    path=Path(path).resolve();cfg=json.loads(path.read_text());require(cfg.get('schema_version')==1,'Unsupported configuration schema');cfg['_file']=str(path);cfg['_base']=str(path.parent)
    if threads is not None:cfg['threads']=threads
    require(isinstance(cfg.get('threads'),int) and cfg['threads']>=1,'threads must be positive')
    require(isinstance(cfg.get('chunk_reads'),int) and cfg['chunk_reads']>0,'chunk_reads must be positive')
    cfg['_out']=str(resolve(cfg,cfg['output']));return cfg

def resolve(cfg,p):
    p=Path(os.path.expandvars(os.path.expanduser(str(p))));return p.resolve() if p.is_absolute() else (Path(cfg['_base'])/p).resolve()

def manifest(cfg,check_files=True):
    rows=read_tsv(resolve(cfg,cfg['samples']));require(rows,'Empty sample manifest');needed={'sample_id','marker','run_id','barcode','fastq'};require(needed<=set(rows[0]),f'Manifest requires: {sorted(needed)}')
    seen={};expanded=[];out=Path(cfg['_out'])
    for row in rows:
        for k in ['sample_id','marker','run_id']:safe_id(row[k])
        require(row['marker'] in cfg['markers'],f'Unknown marker {row["marker"]}')
        require(row['barcode'].lower()!='unclassified','unclassified is not a known biological sample; supply a separate explicit sample ID and provenance if its use is justified')
        pattern=str(resolve(cfg,row['fastq']));files=sorted(Path(f).resolve() for f in glob.glob(pattern,recursive=True)) if check_files else []
        if check_files:require(files,f'No files match: {pattern}')
        for f in files:
            require(f.is_file() and any(str(f).endswith(x) for x in ['.fastq','.fq','.fastq.gz','.fq.gz']),f'Expected FASTQ: {f}')
            require(not f.is_relative_to(out),'Output cannot be a raw input')
            require(f not in seen,f'File repeated in manifest ({seen.get(f)} and {row["sample_id"]}): {f}');seen[f]=row['sample_id']
        expanded.append(dict(row,files=files))
    return expanded

def validate(cfg,check_files=True):
    rows=manifest(cfg,check_files);meta=read_tsv(resolve(cfg,cfg['metadata']));require(meta and 'sample_id' in meta[0],'Metadata needs sample_id');ids=[r['sample_id'] for r in meta];require(len(ids)==len(set(ids)),'Duplicate biological sample IDs in metadata')
    require({r['sample_id'] for r in rows}<=set(ids),'Manifest sample missing in metadata');warnings=[]
    for name,m in cfg['markers'].items():
        require(m['marker'] in {'ssu','its','lsu','long','16s'},f'Unsupported marker {m["marker"]}')
        lo,hi=m.get('min_length'),m.get('max_length');require(isinstance(lo,int) and isinstance(hi,int) and 0<lo<=hi,f'{name}: explicitly select min_length/max_length after checking the amplicon; generic marker profiles have no universal window')
        for v in ['identity','coverage','score_ratio']:require(0<m[v]<=1,f'{name}: invalid {v}')
        require(m.get('trim') in {'none','linked'},f'{name}: trim must be none or linked')
        if m['trim']=='linked':require(m.get('forward') and m.get('reverse'),f'{name}: verified forward/reverse primer sequences required')
        else:warnings.append(f'{name}: no primer trimming; record whether already trimmed or still unverified')
    tasks=json.loads(resolve(cfg,cfg['tasks']).read_text());seen_tasks=set();designs=[]
    for t in tasks:
        safe_id(t['id']);require(t['id'] not in seen_tasks,'Duplicate task ID');seen_tasks.add(t['id'])
        required=[t['group']]+t.get('covariates',[])+([t['strata']] if t.get('strata') else [])+list(t.get('subset',{}));require(set(required)<=set(meta[0]),f'Missing metadata columns for task {t["id"]}')
        require(t['group'] not in t.get('covariates',[]),f'{t["id"]}: group cannot also be a covariate')
        selected=[r for r in meta if all(r.get(k) in (v if isinstance(v,list) else [v]) for k,v in t.get('subset',{}).items())]
        require(selected,f'Empty task subset {t["id"]}')
        formula_columns=[t['group']]+t.get('covariates',[])+([t['strata']] if t.get('strata') else [])
        require(all(re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',col) for col in formula_columns),f'{t["id"]}: formula/stratum column names must use letters, digits and underscores')
        require(all(r.get(col,'')!='' for r in selected for col in formula_columns),f'{t["id"]}: missing group/covariate/stratum values in selected metadata')
        groups={g:sum(r[t['group']]==g for r in selected) for g in sorted({r[t['group']] for r in selected})};designs.append(dict(task=t['id'],N=len(selected),groups=groups))
        if len(groups)<2 or min(groups.values())<2:warnings.append(t['id']+': no replicated group contrast; descriptive distances remain available')

    return {'samples':len({r['sample_id'] for r in rows}),'libraries':len(rows),'raw_files':sum(len(r['files']) for r in rows),'tasks':len(tasks),'task_designs':designs,'warnings':warnings}

def doctor():
    return {name:shutil.which(name) for name in ['python3','minimap2','cutadapt','7z','Rscript','lefse_format_input.py','lefse_run.py']}
