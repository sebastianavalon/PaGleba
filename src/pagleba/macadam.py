import collections,json,re,sqlite3
from pathlib import Path
from .common import require,sha,atomic_json,write_tsv
from .ncbi import Taxonomy
from .functions import ANNOTATION_FIELDS

def import_macadam(database,ncbi,rules,out):
    rules_path=Path(rules);c=sqlite3.connect('file:'+str(Path(database).resolve())+'?mode=ro',uri=True);c.execute('PRAGMA query_only=ON');require(c.execute('PRAGMA quick_check').fetchone()[0]=='ok','MACADAM SQLite failed integrity check');tx=Taxonomy(ncbi);rules=json.loads(Path(rules).read_text());traits=rules['traits'];gates=rules
    tables={r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table'")};require({'taxonomy','faprotax','IJSEMPhenoDB','pathway'}<=tables,'Unsupported MACADAM schema')
    organisms={};bytax=collections.defaultdict(set);names={};rank_byid={};genus_bygenome={}
    for pgdb,path in c.execute('SELECT DISTINCT ID,taxonomy FROM pathway'):
        prefix=pgdb.split('.')[0];ids=([int(prefix)] if prefix.isdigit() else [])+[int(v) for v in reversed(path.strip('.').split('.')) if v.isdigit()];chosen=None
        for tid in ids:
            lin=tx.lineage(tid)
            if any(t in {2,2157} for t,rank,name in lin):chosen=lin;break
        if not chosen:continue
        genus=next((t for t,r,n in reversed(chosen) if r=='genus'),None)
        if not genus:continue
        organisms[pgdb]=chosen;genus_bygenome[pgdb]=genus
        for tid,rank,name in chosen:bytax[tid].add(pgdb);names[tid]=name;rank_byid[tid]=rank
    annotations=[];evidence=[]
    for tid,trait in c.execute('SELECT DISTINCT taxID,pathwayName FROM faprotax'):
        lin=tx.lineage(int(tid))
        if not lin:continue
        tid,rank,name=lin[-1]
        if rank.title() not in ['Domain','Kingdom','Phylum','Class','Order','Family','Genus','Species']:continue
        annotations.append(dict(Source='MACADAM_FAPROTAX',Rank=rank.title(),Taxon=name,Trait=trait,Weight=1,Evidence=f'Explicit reference clade annotation; current taxID {tid}',ContextFamily=''))
    for tid,strain,trait in c.execute("SELECT DISTINCT taxID,strainName,pathway FROM IJSEMPhenoDB WHERE pathway!=''"):
        lin=tx.lineage(int(tid));genus=next(((t,n) for t,r,n in reversed(lin) if r=='genus'),None)
        if not genus:continue
        annotations.append(dict(Source='MACADAM_IJSEM',Rank='Genus',Taxon=genus[1],Trait=trait,Weight=1,Evidence=f'Positive phenotype in reference strain {strain}; extrapolation to genus, not universal',ContextFamily=''))
    lookup=collections.defaultdict(set)
    for trait,rule in traits.items():
        for pid in rule['pathway_ids']:lookup[pid].add(trait)
    positives=collections.defaultdict(set)
    # Score gate and deduplication operate on genome/pathway evidence, never on repeated hierarchy rows.
    seen=set()
    for pgdb,pid,score in c.execute('SELECT ID,pathwayFrameID,pathwayScore FROM pathway WHERE pathwayScore>=?',(rules.get('pathway_score_min',.999999),)):
        if pgdb not in organisms or pid not in lookup or (pgdb,pid) in seen:continue
        seen.add((pgdb,pid))
        for trait in lookup[pid]:positives[trait].add(pgdb)
    for tid,genomes in bytax.items():
        rank=rank_byid[tid]
        if rank.title() not in ['Genus','Family','Order','Class','Phylum','Kingdom','Domain']:continue
        genera={genus_bygenome[g] for g in genomes}
        for trait,pos in positives.items():
            hit=len(genomes&pos);fraction=hit/len(genomes)
            allowed=hit>0 and (fraction>=gates.get('genus_positive_reference_fraction_min',1/3) if rank=='genus' else fraction>=gates.get('higher_rank_positive_reference_fraction_min',.9) and len(genomes)>=gates.get('higher_rank_min_genomes',5) and len(genera)>=gates.get('higher_rank_min_genera',3))
            if hit:evidence.append(dict(TaxID=tid,Rank=rank,Taxon=names[tid],Trait=trait,Reference_Genomes=len(genomes),Positive_Genomes=hit,Reference_Genera=len(genera),Fraction=fraction,Accepted=allowed))
            if allowed:annotations.append(dict(Source='MACADAM_MetaCyc',Rank=rank.title(),Taxon=names[tid],Trait=trait,Weight=1,Evidence=f'Complete pathway; {hit}/{len(genomes)} reference genomes, {len(genera)} genera; gate in supplied rules',ContextFamily=''))
    unique={(r['Source'],r['Rank'],r['Taxon'],r['Trait']):r for r in annotations};write_tsv(out,ANNOTATION_FIELDS,unique.values());write_tsv(str(out)+'.evidence.tsv',['TaxID','Rank','Taxon','Trait','Reference_Genomes','Positive_Genomes','Reference_Genera','Fraction','Accepted'],evidence);atomic_json(str(out)+'.source.json',{'MACADAM_sha256':sha(database),'NCBI_sha256':sha(ncbi),'rules_sha256':sha(rules_path),'reference_genomes':len(organisms),'associations':len(unique),'inference':'reference-based functional potential; absent annotation is unknown; genus splitting resolved via current taxID lineage, not wholesale renaming'})
