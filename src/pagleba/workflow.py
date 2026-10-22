import json,os,subprocess,sys,platform
from pathlib import Path
from . import __version__
from .common import atomic_json,lock,require,read_tsv,sha,write_tsv,RANKS,digest,checked_checkpoint,commit,dependency_versions,r_versions
from .config import load,validate,manifest,resolve,doctor
from .preprocess import merge,trim
from .classify import classify
from .tables import export_tables
from .functions import run_functions
from .statistics import run_statistics
from .lefse import run_lefse
from .report import build

STAGES=['merge','trim','classify','tables','functions','statistics','models','lefse','report']
def existing_tables(cfg):
    if cfg.get('input_tables'):return [(name,resolve(cfg,v['path'])) for name,v in cfg['input_tables'].items() if v['kind']=='taxonomy']
    return [(m,Path(cfg['_out'])/'tables'/m/'analysis_counts.tsv') for m in cfg['markers'] if (Path(cfg['_out'])/'tables'/m/'analysis_counts.tsv').exists()]
def existing_functions(cfg):
    imported=[(name,resolve(cfg,v['path'])) for name,v in cfg.get('input_tables',{}).items() if v['kind']=='traits']
    return imported+[(s['id'],Path(cfg['_out'])/'functions'/s['id']/'trait_fraction.tsv') for s in cfg.get('functions',[]) if (Path(cfg['_out'])/'functions'/s['id']/'trait_fraction.tsv').exists()]

def model_settings(cfg):
    root=Path(cfg['_out']);settings=dict(cfg.get('models',{}));settings.update(jobs=str(root/'statistics/jobs.json'),out=str(root/'models'),permutations=cfg['statistics'].get('permutations',9999),network_bootstrap=cfg['statistics'].get('network_bootstrap',1000),network_rho=cfg['statistics'].get('network_rho',.6),network_q=cfg['statistics'].get('network_q',.05),bootstrap=cfg['statistics'].get('bootstrap',1000),seed=cfg.get('seed',20261021),threads=cfg['threads']);return settings

def run(cfg,stages=None,pilot=None,trees=False):
    imported=bool(cfg.get('input_tables'));stages=stages or ([s for s in STAGES if s not in {'merge','trim','classify','tables'}] if imported else STAGES)
    require(set(stages)<=set(STAGES),'Unknown stage requested');root=Path(cfg['_out']);root.mkdir(parents=True,exist_ok=True)
    if not imported:validate(cfg,check_files='merge' in stages)
    else:
        meta=read_tsv(resolve(cfg,cfg['metadata']));require(meta and 'sample_id' in meta[0],'Metadata needs sample_id');require(len({r['sample_id'] for r in meta})==len(meta),'Duplicate metadata IDs')
        for v in cfg['input_tables'].values():require(resolve(cfg,v['path']).is_file(),'Missing input table')
    with lock(root):
        statepath=root/'workflow.json';state=json.loads(statepath.read_text()) if statepath.exists() else {'stages':{}}
        runtime={'PaGleba':__version__,'dependencies':dependency_versions(['numpy','scipy','matplotlib']),'python':sys.version,'platform':platform.platform(),'commands':doctor(),'config_file':cfg['_file'],'config_sha256':sha(cfg['_file'])};atomic_json(root/'runtime.json',runtime);atomic_json(root/'effective_config.json',{k:v for k,v in cfg.items() if not k.startswith('_')})
        tables=existing_tables(cfg);functions=existing_functions(cfg)
        for stage in stages:
            state['status']='running';state['active_stage']=stage;atomic_json(statepath,state);print('STAGE',stage,flush=True)
            try:
                if not imported and stage in {'functions','statistics','models','lefse','report'}:
                    expected=sorted(c.casefold() for c in cfg.get('taxonomy_filter',{}).get('exclude_clades',[]))
                    for name,path in tables:
                        qc=json.loads((path.parent/'count_reconciliation.json').read_text())
                        require(qc['excluded_clades']==expected,'Configured taxon filter differs from existing analysis tables. Use a new output analysis and rebuild tables; do not relabel old tables in a report.')
                if stage=='merge':merge(cfg,manifest(cfg))
                elif stage=='trim':trim(cfg,json.loads((root/'merge/manifest.json').read_text()))
                elif stage=='classify':
                    result=classify(cfg,json.loads((root/'trim/manifest.json').read_text()),pilot)
                    if result is None:state.update(status='pilot_complete',active_stage=None);atomic_json(statepath,state);return
                elif stage=='tables':tables=export_tables(cfg,json.loads((root/'taxonomy/manifest.json').read_text()))
                elif stage=='functions':functions=run_functions(cfg,tables)+[(n,resolve(cfg,v['path'])) for n,v in cfg.get('input_tables',{}).items() if v['kind']=='traits']
                elif stage=='statistics':run_statistics(cfg,tables,functions)
                elif stage=='models':
                    require((root/'statistics/jobs.json').exists(),'Run statistics first');settings=model_settings(cfg);atomic_json(root/'model_settings.json',settings);folder=root/'models';folder.mkdir(exist_ok=True)
                    signature=digest({'algorithm':'models-v1','R_runtime':r_versions(['jsonlite','igraph','vegan','permute','plspm'],settings.get('r_library')),'settings':settings,'code':sha(Path(__file__).with_name('models.R')),'inputs':{str(p):sha(p) for p in (root/'statistics').rglob('*.tsv')}})
                    if not checked_checkpoint(folder,signature):
                        with (folder/'models.log').open('w') as log:subprocess.run(['Rscript',str(Path(__file__).with_name('models.R')),'--settings',str(root/'model_settings.json')],stdout=log,stderr=subprocess.STDOUT,check=True)
                        commit(folder,signature,[p for p in folder.rglob('*') if p.is_file() and p.name!='COMPLETE.json'])
                elif stage=='lefse':run_lefse(cfg,tables)
                elif stage=='report':build(cfg,tables,functions,trees)
                state['stages'][stage]={'status':'complete'};atomic_json(statepath,state)
            except BaseException as e:
                state['stages'][stage]={'status':'interrupted' if isinstance(e,KeyboardInterrupt) else 'failed','error':str(e)};state['status']='incomplete';atomic_json(statepath,state);raise
        state['status']='requested_stages_complete';state['active_stage']=None;atomic_json(statepath,state);print('COMPLETE:',root)
