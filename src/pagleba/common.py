from __future__ import annotations
import csv,gzip,hashlib,json,os,re,subprocess,time,fcntl
from pathlib import Path
from contextlib import contextmanager

RANKS=['Domain','Kingdom','Phylum','Class','Order','Family','Genus','Species']
class WorkflowError(RuntimeError): pass

def require(test,message):
    if not test: raise WorkflowError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,ensure_ascii=False).encode()).hexdigest()
def atomic_json(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_name(path.name+'.partial');tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n');tmp.replace(path)
def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(path,fields,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_name(path.name+'.partial')
    with tmp.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
    tmp.replace(path)
def gzopen(path,mode='rt'):return gzip.open(path,mode) if str(path).endswith('.gz') else open(path,mode)
def fastq(path):
    with gzopen(path) as f:
        n=0
        while True:
            h=f.readline()
            if not h:break
            s=f.readline().rstrip('\r\n');p=f.readline();q=f.readline().rstrip('\r\n');n+=1
            require(h.startswith('@') and p.startswith('+') and s and len(s)==len(q),f'Invalid four-line FASTQ at record {n}: {path}')
            require(set(s.upper())<=set('ACGTRYSWKMBDHVNU'),f'Invalid nucleotide sequence at record {n}: {path}')
            require(all(33<=ord(c)<=126 for c in q),f'Invalid FASTQ quality at record {n}: {path}')
            yield h.rstrip('\r\n'),s,p.rstrip('\r\n'),q

def safe_id(s):require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',s)) and s not in {'.','..'},f'Unsafe identifier: {s!r}');return s

def command(argv,log=None):
    argv=[str(v) for v in argv]
    if log:
        with Path(log).open('w') as f:subprocess.run(argv,stdout=f,stderr=subprocess.STDOUT,check=True)
    else:subprocess.run(argv,check=True)

def checked_checkpoint(folder,signature):
    folder=Path(folder);marker=folder/'COMPLETE.json'
    if not marker.exists():return False
    state=json.loads(marker.read_text());require(state['signature']==signature,f'Changed input/configuration in completed stage {folder}. Use a new --out/project.')
    for name,value in state['files'].items():require((folder/name).is_file() and sha(folder/name)==value,f'Corrupt completed output: {folder/name}')
    return True

def commit(folder,signature,files,details=None):
    folder=Path(folder);atomic_json(folder/'COMPLETE.json',{'signature':signature,'files':{str(Path(p).relative_to(folder)):sha(p) for p in files},'details':details or {},'completed_epoch':time.time()})

@contextmanager
def lock(folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    with (folder/'runner.lock').open('a') as f:
        try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise WorkflowError(f'Another process owns {folder}')
        try:yield
        finally:fcntl.flock(f,fcntl.LOCK_UN)

def dependency_versions(names):
    import importlib.metadata
    result={}
    for name in names:
        try:result[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:result[name]='not_installed'
    return result

def r_versions(packages,library_path=None):
    import shutil
    if not shutil.which('Rscript'):return 'Rscript_not_installed'
    libs=('.libPaths(c('+json.dumps(str(library_path))+',.libPaths()));') if library_path else ''
    code=libs+'cat("R=",as.character(getRversion()),"\\n",sep="");for(p in c('+','.join(json.dumps(p) for p in packages)+'))cat(p,"=",if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else"not_installed","\\n",sep="")'
    return subprocess.check_output(['Rscript','-e',code],text=True).strip()
