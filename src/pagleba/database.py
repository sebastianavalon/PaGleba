"""Reference preparation without loading the whole EUKARYOME taxonomy into RAM."""
import csv,gzip,json,os,re,shutil,sqlite3,subprocess,urllib.request,zipfile
from pathlib import Path
from .common import RANKS,require,sha,atomic_json,gzopen,command,lock,WorkflowError

def fasta_records(path):
    with gzopen(path) as f:
        name=None;seq=[]
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(seq)
                name=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
        if name is not None:yield name,''.join(seq)

def unpack(archive,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as z:
            require(z.testzip() is None,'Corrupt ZIP')
            for info in z.infolist():
                if info.is_dir():continue
                require(info.file_size<30*1024**3,'Unexpected archive member size')
                target=out/Path(info.filename).name
                require(target.name not in {'','..'},'Unsafe archive member')
                with z.open(info) as a,target.open('wb') as b:shutil.copyfileobj(a,b)
        for arc in out.glob('*.7z'):
            require(shutil.which('7z'),'Install 7z for nested EUKARYOME archive')
            command(['7z','x','-y','-o'+str(out/'nested'),arc])
    elif str(archive).endswith('.7z'):
        require(shutil.which('7z'),'Install 7z');command(['7z','x','-y','-o'+str(out),archive])
    else:raise WorkflowError('Use a ZIP/7z archive or explicit --fasta and --taxonomy')
    return out

def fetch_eukaryome(marker,version,out,url=None,expected=None):
    require(marker in {'ssu','its','lsu','long'},'EUKARYOME marker must be ssu/its/lsu/long')
    require(bool(re.fullmatch(r'\d+\.\d+(?:\.\d+)?',version)),'Explicit numeric release version required')
    token={'ssu':'SSU','its':'ITS','lsu':'LSU','long':'long'}[marker];out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    url=url or f'https://sisu.ut.ee/wp-content/uploads/sites/643/mothur_EUK_{token}_v{version}.zip';dest=out/f'mothur_EUK_{token}_v{version}.zip'
    if not dest.exists():
        tmp=dest.with_suffix('.partial')
        with urllib.request.urlopen(url,timeout=120) as r,tmp.open('wb') as f:shutil.copyfileobj(r,f,8*1024*1024)
        tmp.replace(dest)
    actual=sha(dest);require(not expected or actual==expected,'Publisher/user SHA256 mismatch')
    atomic_json(out/'download.json',{'url':url,'provider':'eukaryome','marker':marker,'version':version,'sha256':actual,'expected_sha256':expected,'license':'CC-BY-4.0; reference data not bundled with code'})
    ex=unpack(dest,out/'unpacked');fa=list(ex.rglob(f'mothur_EUK_{token}_v{version}.fasta'));tax=list(ex.rglob(f'mothur_EUK_{token}_v{version}.tax'));require(len(fa)==len(tax)==1,'Reference archive layout differs: inspect contents and pass explicit FASTA/taxonomy')
    return fa[0],tax[0]

def tax_rows(path,fmt,provider,rank_map=None):
    if fmt=='mothur':
        require(provider=='eukaryome','Positional mothur ranks are only defined here for EUKARYOME, not generic bacterial mothur files')
        with gzopen(path) as f:
            for line in f:
                if not line.strip():continue
                rid,tax=line.rstrip('\r\n').split('\t',1);levels=[re.sub(r'\(\d+(?:\.\d+)?\)$','',v) for v in tax.strip(';').split(';')];require(1<=len(levels)<=7,'Unexpected EUKARYOME rank count: use explicit rank TSV if schema changes')
                domain='Bacteria' if levels[0]=='_Bacteria' else 'Archaea' if levels[0]=='_Archaea' else 'Eukaryota'
                yield rid,([domain]+levels+['']*8)[:8],tax
    elif fmt=='silva':
        require(rank_map,'SILVA requires --rank-map tax_slv_ssu_*.txt(.gz)')
        rank_paths={}
        with gzopen(rank_map) as f:
            for line in f:
                v=line.rstrip('\r\n').split('\t')
                if len(v)>=3 and v[1].isdigit():rank_paths[v[0]]=v[2].casefold()
        with gzopen(path) as f:
            for row in csv.DictReader(f,delimiter='\t'):
                rid='.'.join(row[k] for k in ['primaryAccession','start','stop']);original=row['path'];parts=original.strip(';').split(';');levels=['']*8
                for i,val in enumerate(parts):
                    prefix=';'.join(parts[:i+1])+';';rank=rank_paths.get(prefix,'');rank='domain' if rank in {'superkingdom','domain'} else rank
                    if rank.title() in RANKS:levels[RANKS.index(rank.title())]=val
                require(levels[0],f'Missing SILVA domain for {rid}');yield rid,levels,original
    else:
        with gzopen(path) as f:
            for row in csv.DictReader(f,delimiter='\t'):
                rid=row.get('reference') or row.get('Reference')
                require(rid and any(k in row for k in RANKS),'Taxonomy TSV needs reference and explicit Domain..Species columns')
                levels=[row.get(r,'').strip() for r in RANKS];require(levels[0],f'Missing domain for {rid}');yield rid,levels,';'.join(levels)+';'

def prepare(fasta,taxonomy,out,provider,marker,version,fmt='mothur',rank_map=None,preset='map-ont',threads=8):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True);fasta=Path(fasta).resolve();taxonomy=Path(taxonomy).resolve()
    require(not (out/'database.json').exists(),'Prepared database already exists. Choose a new --out for any change.')
    require(shutil.which('minimap2'),'minimap2 not on PATH');require(marker in {'ssu','its','lsu','long','16s'},'Unknown marker')
    with lock(out):
        sql=out/'taxonomy.sqlite.partial';sql.unlink(missing_ok=True);c=sqlite3.connect(sql);c.execute('PRAGMA journal_mode=DELETE');c.execute('CREATE TABLE refs(reference TEXT PRIMARY KEY, levels TEXT NOT NULL, original TEXT, seen INTEGER DEFAULT 0)')
        for rid,levels,original in tax_rows(taxonomy,fmt,provider,rank_map):
            c.execute('INSERT INTO refs VALUES(?,?,?,0)',(rid,json.dumps(levels),original))
        c.commit();n_tax=c.execute('SELECT count(*) FROM refs').fetchone()[0];require(n_tax>0,'Empty taxonomy')
        count=bases=0;dna=out/'reference.fasta.partial'
        with dna.open('w') as f:
            for rid,seq in fasta_records(fasta):
                seq=seq.upper().replace('U','T').replace('-','').replace('.','');require(seq and set(seq)<=set('ACGTRYSWKMBDHVN'),f'Invalid reference alphabet {rid}')
                cur=c.execute('UPDATE refs SET seen=1 WHERE reference=? AND seen=0',(rid,));require(cur.rowcount==1,f'Missing or duplicate FASTA ID {rid}')
                f.write('>'+rid+'\n'+seq+'\n');count+=1;bases+=len(seq)
                if count%10000==0:c.commit()
        c.commit();missing=c.execute('SELECT count(*) FROM refs WHERE seen=0').fetchone()[0]
        # SILVA taxmap may contain additional reference records; their count is explicit.
        require(provider=='silva' or missing==0,f'{missing} taxonomy IDs absent from FASTA; incompatible input pair')
        c.execute('DELETE FROM refs WHERE seen=0');c.commit();c.close();sql.replace(out/'taxonomy.sqlite');dna.replace(out/'reference.fasta')
        temp=out/'reference.partial.mmi';cmd=['minimap2','-x',preset,'-t',str(threads),'-I',str(max(bases+100000,1000000)),'-d',str(temp),str(out/'reference.fasta')];command(cmd,out/'index.log');temp.replace(out/'reference.mmi')
        meta={'schema_version':1,'provider':provider,'marker':marker,'version':version,'preset':preset,'sequences':count,'bases':bases,'unselected_taxmap_records':missing,'ranks':RANKS,'input_sha256':{'fasta':sha(fasta),'taxonomy':sha(taxonomy)},'files':{n:sha(out/n) for n in ['reference.fasta','reference.mmi','taxonomy.sqlite']},'minimap2_version':subprocess.check_output(['minimap2','--version'],text=True).strip(),'index_command':cmd,'taxonomy_schema':'EUKARYOME positional kingdom..species plus local domain' if fmt=='mothur' else 'explicit ranked taxonomy','taxonomic_names':'database release names retained; not independently updated to latest nomenclature'}
        if rank_map:meta['rank_map_sha256']=sha(rank_map)
        atomic_json(out/'database.json',meta)
    return meta
