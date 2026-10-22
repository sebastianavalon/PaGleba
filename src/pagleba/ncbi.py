"""Local same-taxID nomenclature reconciliation. A genus split is never a synonym substitution."""
import functools,json,sqlite3
from pathlib import Path
from .common import require,sha,atomic_json,read_tsv,write_tsv,RANKS

def prepare(nodes,names,out,merged=None):
    out=Path(out).resolve();require(not out.exists(),f'NCBI SQLite already exists: {out}');out.parent.mkdir(parents=True,exist_ok=True);tmp=out.with_suffix('.partial');tmp.unlink(missing_ok=True);c=sqlite3.connect(tmp)
    c.executescript('CREATE TABLE taxon(taxid INTEGER PRIMARY KEY,parent INTEGER,rank TEXT,name TEXT);CREATE TABLE alias(name TEXT,taxid INTEGER,kind TEXT,original_name TEXT);CREATE TABLE merged(old INTEGER PRIMARY KEY,new INTEGER);')
    batch=[]
    with Path(nodes).open() as f:
        for line in f:
            v=[p.strip() for p in line.split('|')];batch.append((int(v[0]),int(v[1]),v[2]))
            if len(batch)>=10000:c.executemany('INSERT INTO taxon(taxid,parent,rank) VALUES(?,?,?)',batch);batch=[]
    if batch:c.executemany('INSERT INTO taxon(taxid,parent,rank) VALUES(?,?,?)',batch)
    c.commit();batch=[];updates=[]
    with Path(names).open() as f:
        for line in f:
            v=[p.strip() for p in line.split('|')];tid=int(v[0]);name=v[1];kind=v[3]
            if kind=='scientific name':updates.append((name,tid))
            if kind in {'scientific name','synonym','equivalent name','genbank synonym'}:batch.append((name.casefold(),tid,kind,name))
            if len(batch)>=20000:
                c.executemany('INSERT INTO alias VALUES(?,?,?,?)',batch);c.executemany('UPDATE taxon SET name=? WHERE taxid=?',updates);c.commit();batch=[];updates=[]
    c.executemany('INSERT INTO alias VALUES(?,?,?,?)',batch);c.executemany('UPDATE taxon SET name=? WHERE taxid=?',updates)
    if merged:
        with Path(merged).open() as f:
            c.executemany('INSERT INTO merged VALUES(?,?)',[(int(v[0]),int(v[1])) for line in f if (v:=[p.strip() for p in line.split('|')])])
    c.execute('CREATE INDEX alias_name ON alias(name)');c.commit();c.close();tmp.replace(out);atomic_json(str(out)+'.json',{'nodes_sha256':sha(nodes),'names_sha256':sha(names),'merged_sha256':sha(merged) if merged else None,'method':'scientific/synonym/equivalent names within the same stable taxID; merged IDs followed; no fuzzy genus reassignment'})

class Taxonomy:
    def __init__(self,path):self.c=sqlite3.connect('file:'+str(Path(path).resolve())+'?mode=ro',uri=True)
    @functools.lru_cache(maxsize=100000)
    def current(self,tid):
        seen=set()
        while tid not in seen:
            seen.add(tid);r=self.c.execute('SELECT new FROM merged WHERE old=?',(tid,)).fetchone()
            if not r:break
            tid=r[0]
        return tid
    @functools.lru_cache(maxsize=100000)
    def lineage(self,tid):
        tid=self.current(int(tid));out=[];seen=set()
        while tid not in seen:
            seen.add(tid);r=self.c.execute('SELECT parent,rank,name FROM taxon WHERE taxid=?',(tid,)).fetchone()
            if not r:break
            out.append((tid,r[1],r[2] or ''))
            if r[0]==tid:break
            tid=r[0]
        return tuple(reversed(out))
    def resolve(self,name,rank):
        choices={self.current(r[0]) for r in self.c.execute('SELECT taxid FROM alias WHERE name=?',(name.casefold(),))};out=[]
        for tid in choices:
            lin=self.lineage(tid)
            if lin and lin[-1][1] in ({'domain','superkingdom'} if rank.casefold()=='domain' else {rank.casefold()}):out.append(lin[-1])
        return out

def aliases_for_table(table,annotations,ncbi,out):
    tx=Taxonomy(ncbi);names=set()
    for r in read_tsv(table):names.update((rank,r.get(rank,'')) for rank in RANKS if r.get(rank,''))
    for r in read_tsv(annotations):
        if r.get('Rank','').title() in RANKS:names.add((r['Rank'].title(),r['Taxon']))
    rows=[]
    for rank,name in sorted(names):
        got=tx.resolve(name,rank)
        if len(got)==1:
            tid,_,canonical=got[0];rows.append(dict(Rank=rank,Original=name,Canonical=canonical,TaxID=tid,Evidence='NCBI same taxID scientific name/synonym; '+str(tid)))
    write_tsv(out,['Rank','Original','Canonical','TaxID','Evidence'],rows)
