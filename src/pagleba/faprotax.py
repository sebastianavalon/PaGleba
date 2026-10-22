"""Adapter to the user's official FAPROTAX collapse_table.py, not a replacement rules parser."""
import csv,json,subprocess,sys
from pathlib import Path
from .common import RANKS,read_tsv,write_tsv,require,sha,atomic_json
from .functions import usable,ANNOTATION_FIELDS

def import_official(table,reference,script,out,python=None,ncbi=None):
    rows=read_tsv(table);require(rows,'Empty FAPROTAX input');folder=Path(str(out)+'.work');folder.mkdir(parents=True,exist_ok=True);taxonomies=[]
    from .ncbi import Taxonomy
    tx=Taxonomy(ncbi) if ncbi else None;alias_audit=[]
    for r in rows:
        tokens=[]
        for rank in RANKS:
            name=r.get(rank,'')
            if not usable(name):continue
            tokens.append(name)
            if tx:
                resolved=tx.resolve(name,rank)
                if len(resolved)==1:
                    tid,_,canonical=resolved[0];synonyms=[v[0] for v in tx.c.execute('SELECT original_name FROM alias WHERE taxid=?',(tid,))];tokens.extend(synonyms);alias_audit.append(dict(Feature=r['Feature'],Rank=rank,Original=name,TaxID=tid,Added_Synonyms=json.dumps(synonyms)))
        taxonomies.append(';'.join(dict.fromkeys(tokens)))
    write_tsv(folder/'name_compatibility.tsv',['Feature','Rank','Original','TaxID','Added_Synonyms'],alias_audit)
    write_tsv(folder/'input.tsv',['Feature','taxonomy','occurrence'],[dict(Feature=r['Feature'],taxonomy=tax,occurrence=1) for r,tax in zip(rows,taxonomies)])
    cmd=[python or sys.executable,str(script),'-i',str(folder/'input.tsv'),'-g',str(reference),'-o',str(folder/'collapsed.tsv'),'-r',str(folder/'report.txt'),'--out_groups2records_table',str(folder/'mapping.tsv'),'-d','taxonomy','--omit_columns','0','-n','none','--force'];atomic_json(folder/'command.json',cmd)
    with (folder/'run.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    with (folder/'mapping.tsv').open() as f:assigned=list(csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t'))
    require(len(assigned)==len(rows),'FAPROTAX mapping row count mismatch');outrows=[]
    for r,a,tax in zip(rows,assigned,taxonomies):
        require(a.get('record')==tax,'FAPROTAX record-order/identity mismatch')
        for trait,val in a.items():
            if trait!='record' and float(val)>0:outrows.append(dict(Source='FAPROTAX_official',Rank='Feature',Taxon=r['Feature'],Trait=trait,Weight=1,Evidence='Official collapse_table.py match to the complete retained taxonomy',ContextFamily=''))
    write_tsv(out,ANNOTATION_FIELDS,outrows);atomic_json(str(out)+'.source.json',{'FAPROTAX_sha256':sha(reference),'NCBI_sha256':sha(ncbi) if ncbi else None,'official_script_sha256':sha(script),'input_table_sha256':sha(table),'functional_associations':len(outrows)})
