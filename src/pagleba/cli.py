import argparse,collections,glob,json,os,re,sys,traceback,subprocess
from pathlib import Path
from . import __version__
from .common import WorkflowError,require,read_tsv,write_tsv,atomic_json,RANKS
from . import config,database,workflow

def configuration_flag(parser):parser.add_argument('--config',required=True,type=Path,help='Project JSON; relative data paths resolve from this file');parser.add_argument('--threads',type=int,help='Override CPU threads without editing code')
def parser():
    p=argparse.ArgumentParser(prog='pagleba',description='Portable ONT 16S/18S/ITS amplicon taxonomy and ecological analysis. No implicit plant removal.');p.add_argument('--version',action='version',version=__version__);s=p.add_subparsers(dest='cmd',required=True)
    q=s.add_parser('demo',help='Create a tiny synthetic two-marker project, including reference indexes');q.add_argument('--project',required=True,type=Path);q.add_argument('--samples',type=int,default=12);q.add_argument('--threads',type=int,default=2)
    q=s.add_parser('init',help='Create project configuration and manifest templates');q.add_argument('--project',required=True,type=Path);q.add_argument('--profile',choices=list(config.PROFILES)+['pagleba'],default='18s');q.add_argument('--interactive',action='store_true',help='Ask simple setup questions in the terminal')
    q=s.add_parser('set',help='Set a typed JSON configuration value using a dotted key');q.add_argument('--config',required=True,type=Path);q.add_argument('--key',required=True);q.add_argument('--value',required=True,help='JSON literal: 150, 0.95, true, [], or a quoted string')
    for name in ['validate','plan','status']:
        q=s.add_parser(name);configuration_flag(q)
        if name=='validate':q.add_argument('--no-files',action='store_true',help='Check schema/design without expanding input globs')
    q=s.add_parser('doctor',help='Inspect tool availability and optional R dependencies');q.add_argument('--r-check',action='store_true')
    for name in ['run','merge','trim','classify','tables','functions','statistics','models','lefse','report']:
        q=s.add_parser(name,help='Run '+('the configured pipeline' if name=='run' else name+' stage')+'; completed read checkpoints are verified and reused');configuration_flag(q)
        if name=='run':q.add_argument('--stages',help='Comma-separated ordered stages; prerequisites must exist');q.add_argument('--dry-run',action='store_true')
        q.add_argument('--pilot-chunks',type=int,help='Stop after N newly classified chunks; no final table for a partial run');q.add_argument('--trees',action='store_true',help='Include paginated vector family taxonomic hierarchies in report')
    q=s.add_parser('metadata',help='Create a metadata template from an explicit sample manifest');q.add_argument('--samples',required=True,type=Path);q.add_argument('--out',required=True,type=Path)
    q=s.add_parser('samples',help='Build an explicit manifest from barcode directories');ss=q.add_subparsers(dest='samplecmd',required=True);q=ss.add_parser('discover');q.add_argument('--run',action='append',required=True,help='Repeat RUN_ID=/arbitrary/root; roots with the same barcode must represent the same biological sample');q.add_argument('--marker',required=True);q.add_argument('--barcodes',required=True,help='e.g. 1-24 or 1-8,11,15-19');q.add_argument('--sample-prefix',default='');q.add_argument('--rename-start',type=int,help='For a different sample set, assign consecutive new barcode numbers; headers updated on merge');q.add_argument('--out',required=True,type=Path)
    q=s.add_parser('db',help='Prepare marker-specific reference FASTA, SQLite taxonomy and minimap2 index');ds=q.add_subparsers(dest='dbcmd',required=True)
    q=ds.add_parser('eukaryome',help='Download pinned EUKARYOME archive and prepare its matching marker');q.add_argument('--marker',choices=['18s','ssu','its','lsu','long'],required=True);q.add_argument('--version',default='2.1');q.add_argument('--out',required=True,type=Path);q.add_argument('--url');q.add_argument('--sha256',help='Expected trusted archive SHA256 if available');q.add_argument('--preset',default='map-ont');q.add_argument('--threads',type=int,default=8)
    q=ds.add_parser('prepare',help='Prepare local FASTA/taxonomy, including SILVA');q.add_argument('--fasta',required=True,type=Path);q.add_argument('--taxonomy',required=True,type=Path);q.add_argument('--provider',choices=['eukaryome','silva','custom'],required=True);q.add_argument('--marker',choices=['18s','ssu','its','lsu','long','16s'],required=True);q.add_argument('--version',required=True);q.add_argument('--format',choices=['mothur','tsv','silva'],default='mothur');q.add_argument('--rank-map',type=Path);q.add_argument('--out',required=True,type=Path);q.add_argument('--preset',default='map-ont');q.add_argument('--threads',type=int,default=8)
    q=s.add_parser('traits',help='Prepare or reconcile functional reference annotations');ts=q.add_subparsers(dest='traitcmd',required=True)
    q=ts.add_parser('import');q.add_argument('--format',choices=['long','eukfunc','fungaltraits'],required=True);q.add_argument('--reference',required=True,type=Path);q.add_argument('--out',required=True,type=Path);q.add_argument('--fields',help='Comma-separated FungalTraits fields')
    q=ts.add_parser('faprotax');q.add_argument('--table',required=True,type=Path);q.add_argument('--reference',required=True,type=Path);q.add_argument('--script',required=True,type=Path);q.add_argument('--python',help='Interpreter compatible with the official FAPROTAX script');q.add_argument('--ncbi',type=Path,help='Optional same-taxID scientific/synonym bridge for historical FAPROTAX names');q.add_argument('--out',required=True,type=Path)
    q=ts.add_parser('macadam');q.add_argument('--database',required=True,type=Path);q.add_argument('--ncbi',required=True,type=Path);q.add_argument('--rules',required=True,type=Path);q.add_argument('--out',required=True,type=Path)
    q=ts.add_parser('aliases');q.add_argument('--table',required=True,type=Path);q.add_argument('--annotations',required=True,type=Path);q.add_argument('--ncbi',required=True,type=Path);q.add_argument('--out',required=True,type=Path)
    q=s.add_parser('ncbi',help='Build an offline NCBI taxonomy/name reconciliation SQLite');q.add_argument('--nodes',required=True,type=Path);q.add_argument('--names',required=True,type=Path);q.add_argument('--merged',type=Path);q.add_argument('--out',required=True,type=Path)
    q=s.add_parser('anchors',help='Compute explicit linked anchor/target distances, such as the PaGleba B4 question');q.add_argument('--table',action='append',default=[]);q.add_argument('--traits',action='append',default=[]);q.add_argument('--units',required=True,type=Path);q.add_argument('--out',required=True,type=Path);q.add_argument('--relative-threshold',type=float,default=.0001);q.add_argument('--prevalence',type=float,default=.1)
    q=s.add_parser('analyze',help='Start at existing taxonomy/trait tables, bypassing FASTQ stages');q.add_argument('--table',action='append',default=[],help='Repeat BLOCK=/path/counts.tsv (Feature,Domain..Species,samples)');q.add_argument('--traits',action='append',default=[],help='Repeat BLOCK=/path/fraction.tsv (Feature,Source,Trait,samples)');q.add_argument('--metadata',required=True,type=Path);q.add_argument('--tasks',required=True,type=Path);q.add_argument('--out',required=True,type=Path);q.add_argument('--threads',type=int,default=12);q.add_argument('--permutations',type=int,default=9999);q.add_argument('--models-config',type=Path);q.add_argument('--observed-counts',action='store_true',help='Assert taxonomy inputs are unscaled observed integer reads; enables alpha diversity and rarefaction');q.add_argument('--trees',action='store_true');q.add_argument('--stages',default='statistics,models,report')
    return p

def main(argv=None):
    a=parser().parse_args(argv)
    try:
        if a.cmd=='demo':
            from .demo import create
            print(create(a.project,a.samples,a.threads))
        elif a.cmd=='init':
            profile=a.profile
            if a.interactive:
                choice=input('Marker profile [18s / its / 16s / pagleba] (default '+profile+'): ').strip();profile=choice or profile;require(profile in list(config.PROFILES)+['pagleba'],'Unknown profile')
            f=config.initialize(a.project,'18s-f566-r1200' if profile=='pagleba' else profile)
            if profile=='pagleba':
                c=json.loads(f.read_text());m=dict(config.PROFILES['16s']);m.update(database='databases/16s',trim='none',forward=None,reverse=None,error_rate=.15,overlap=16,preset='map-ont',secondary=50);c['markers']['16s']=m;atomic_json(f,c);rows=read_tsv(f.parent/'samples.tsv');rows.append(dict(rows[0],marker='16s',fastq='/replace/with/bacterial_run/barcode01/*.fastq.gz'));write_tsv(f.parent/'samples.tsv',list(rows[0]),rows);tasks=json.loads((f.parent/'tasks.json').read_text());tasks[0]['blocks']=['16s','ssu'];atomic_json(f.parent/'tasks.json',tasks)
            print(f)
        elif a.cmd=='set':
            c=json.loads(a.config.read_text());keys=a.key.split('.');target=c
            for k in keys[:-1]:target=target.setdefault(k,{})
            try:value=json.loads(a.value)
            except json.JSONDecodeError:value=a.value
            target[keys[-1]]=value;atomic_json(a.config,c);print(a.key,'=',value)
        elif a.cmd=='doctor':
            result=config.doctor()
            if a.r_check and result['Rscript']:
                z=subprocess.run(['Rscript','-e','for(p in c("vegan","permute","jsonlite","igraph","plspm"))cat(p,requireNamespace(p,quietly=TRUE),"\\n")'],capture_output=True,text=True);result['R_packages']=z.stdout;result['R_errors']=z.stderr
            print(json.dumps(result,indent=2))
        elif a.cmd=='metadata':
            rows=read_tsv(a.samples);ids=sorted({r['sample_id'] for r in rows});require(not a.out.exists(),'Output exists; choose another metadata file');write_tsv(a.out,['sample_id','group','location','material','stage'],[dict(sample_id=s,group='',location='',material='',stage='') for s in ids]);print(a.out)
        elif a.cmd=='samples':
            wanted=set()
            for part in a.barcodes.split(','):
                v=part.split('-');wanted.update(range(int(v[0]),int(v[-1])+1))
            require(all(n>0 for n in wanted),'Barcodes must be positive');mapping={n:(a.rename_start+i if a.rename_start is not None else n) for i,n in enumerate(sorted(wanted))};rows=[]
            for run in a.run:
                rid,root=run.split('=',1);root=Path(root).expanduser().resolve();require(root.is_dir(),f'Missing run root {root}')
                for f in sorted(root.rglob('*')):
                    if not f.is_file() or not any(str(f).endswith(e) for e in ['.fastq','.fq','.fastq.gz','.fq.gz']):continue
                    matches=[re.fullmatch(r'barcode(\d+)',p.name,re.I) for p in f.parents if p!=root.parent];numbers={int(m.group(1)) for m in matches if m}
                    require(len(numbers)<=1,f'Ambiguous nested barcode directories: {f}')
                    if not numbers or next(iter(numbers)) not in wanted:continue
                    n=next(iter(numbers));source=f'barcode{n:02d}';target=f'barcode{mapping[n]:02d}';rows.append(dict(sample_id=a.sample_prefix+target,marker=a.marker,run_id=rid,barcode=source,fastq=str(f),header_barcode=target))
            require(rows,'No selected barcode FASTQs found');require(not a.out.exists(),'Manifest already exists; choose another --out');write_tsv(a.out,['sample_id','marker','run_id','barcode','fastq','header_barcode'],rows);print('Manifest rows:',len(rows))
        elif a.cmd=='db':
            marker='ssu' if a.marker=='18s' else a.marker
            if a.dbcmd=='eukaryome':fa,tax=database.fetch_eukaryome(marker,a.version,a.out/'download',a.url,a.sha256);result=database.prepare(fa,tax,a.out,'eukaryome',marker,a.version,preset=a.preset,threads=a.threads)
            else:result=database.prepare(a.fasta,a.taxonomy,a.out,a.provider,marker,a.version,a.format,a.rank_map,a.preset,a.threads)
            print(json.dumps(result,indent=2))
        elif a.cmd=='traits':
            from .functions import import_traits
            from .macadam import import_macadam
            from .ncbi import aliases_for_table
            if a.traitcmd=='import':import_traits(a.format,a.reference,a.out,a.fields.split(',') if a.fields else None)
            elif a.traitcmd=='faprotax':
                from .faprotax import import_official
                import_official(a.table,a.reference,a.script,a.out,a.python,a.ncbi)
            elif a.traitcmd=='macadam':import_macadam(a.database,a.ncbi,a.rules,a.out)
            else:aliases_for_table(a.table,a.annotations,a.ncbi,a.out)
        elif a.cmd=='ncbi':
            from .ncbi import prepare
            prepare(a.nodes,a.names,a.out,a.merged)
        elif a.cmd=='anchors':
            from .anchors import run
            parse=lambda values:[(v.split('=',1)[0],Path(v.split('=',1)[1]).resolve()) for v in values]
            require(a.table or a.traits,'Provide at least one block');run(parse(a.table),parse(a.traits),a.units,a.out,a.relative_threshold,a.prevalence)
        elif a.cmd=='analyze':
            require(a.table or a.traits,'Provide at least one input block');root=a.out.resolve();root.mkdir(parents=True,exist_ok=True);f=root/'analysis_config.json';require(not f.exists(),'Analysis config exists. Resume with pagleba run --config analysis_config.json; do not replace it.')
            blocks={}
            for kind,values in [('taxonomy',a.table),('traits',a.traits)]:
                for v in values:
                    name,path=v.split('=',1);require(name not in blocks,'Duplicate block name');blocks[name]={'kind':kind,'path':str(Path(path).expanduser().resolve())}
            c={'schema_version':1,'project':'table_analysis','output':'.','metadata':str(a.metadata.resolve()),'tasks':str(a.tasks.resolve()),'threads':a.threads,'chunk_reads':5000,'seed':20261021,'markers':{},'input_tables':blocks,'counts_provenance':'observed' if a.observed_counts else 'unknown','functions':[],'taxonomy_filter':{'exclude_clades':[]},'statistics':{'permutations':a.permutations,'relative_threshold':.0001,'prevalence':.1,'bootstrap':1000,'network_rho':.6,'network_q':.05}}
            if a.models_config:c['models']=json.loads(a.models_config.read_text())
            atomic_json(f,c);workflow.run(config.load(f),a.stages.split(','),trees=a.trees)
        else:
            cfg=config.load(a.config,a.threads)
            if a.cmd=='validate':print(json.dumps(config.validate(cfg,not a.no_files),indent=2))
            elif a.cmd=='status':
                p=Path(cfg['_out'])/'workflow.json';print(p.read_text() if p.exists() else 'No run state yet')
            elif a.cmd=='plan' or (a.cmd=='run' and a.dry_run):
                if not cfg.get('input_tables'):print(json.dumps(config.validate(cfg),indent=2))
                print(json.dumps({'output':cfg['_out'],'markers':cfg.get('markers'),'input_tables':cfg.get('input_tables'),'stages':workflow.STAGES,'taxonomy_filter':cfg.get('taxonomy_filter'),'note':'No data are modified by plan. Stage-specific external tools are required only when that stage runs.'},indent=2))
            else:
                stages=(a.stages.split(',') if a.stages else None) if a.cmd=='run' else [a.cmd];workflow.run(cfg,stages,a.pilot_chunks,a.trees)
    except KeyboardInterrupt:print('Interrupted. Completed checkpoints are preserved. Run the same command to resume.',file=sys.stderr);return 130
    except subprocess.CalledProcessError as e:
        print('ERROR: external tool failed. Inspect the active stage log in the configured output directory. Command:',e.cmd,file=sys.stderr);return e.returncode or 2
    except (WorkflowError,FileNotFoundError,ValueError,KeyError) as e:print('ERROR:',e,file=sys.stderr);return 2
    return 0
if __name__=='__main__':sys.exit(main())
