import collections,gzip,json,re,shutil,subprocess
from pathlib import Path
from .common import sha,digest,checked_checkpoint,commit,fastq,require,atomic_json,command,write_tsv

BARCODE=re.compile(r'(?<![A-Za-z0-9])barcode\d+(?![A-Za-z0-9])',re.I)
def rewrite_description(line,source,target):
    if source==target or not target:return line,0
    first,sep,rest=line.partition(' ')
    if not sep:return line,0
    n=0
    def change(m):
        nonlocal n
        if m.group().casefold()==source.casefold():n+=1;return target
        return m.group()
    return first+sep+BARCODE.sub(change,rest),n

def merge(cfg,rows):
    out=Path(cfg['_out'])/'merge';grouped=collections.defaultdict(list)
    for row in rows:grouped[(row['marker'],row['sample_id'])].append(row)
    result=[]
    for (marker,sample),entries in sorted(grouped.items()):
        folder=out/marker/sample;folder.mkdir(parents=True,exist_ok=True)
        sig=digest({'algorithm':'merge-v1','code_sha256':sha(Path(__file__)),'inputs':[{**{k:v for k,v in r.items() if k!='files'},'files':{str(f):sha(f) for f in r['files']}} for r in entries]})
        dest=folder/(sample+'.fastq.gz')
        if checked_checkpoint(folder,sig):print('RESUME merge',marker,sample,flush=True)
        else:
            temporary=folder/(sample+'.partial.fastq.gz');reads=changed=missing=0
            with gzip.open(temporary,'wt',compresslevel=1) as f:
                for row in entries:
                    source=row['barcode'];target=row.get('header_barcode') or source
                    for sourcefile in row['files']:
                        for h,s,p,q in fastq(sourcefile):
                            h,n=rewrite_description(h,source,target);p,_=rewrite_description(p,source,target)
                            if source!=target and not n:missing+=1
                            changed+=bool(n);reads+=1;f.write(h+'\n'+s+'\n'+p+'\n'+q+'\n')
            require(not missing or cfg.get('allow_missing_header_barcode',False),f'{sample}: {missing} renamed reads lack a matching barcode in the header description. Check manifest or explicitly set allow_missing_header_barcode=true; UUIDs will never be changed.')
            require(reads>0,f'Empty merged sample {sample}');temporary.replace(dest)
            atomic_json(folder/'qc.json',{'reads':reads,'header_records_renamed':changed,'header_records_without_matching_barcode':missing,'duplicates':'all occurrences retained','source_files':[str(f) for r in entries for f in r['files']]})
            commit(folder,sig,[dest,folder/'qc.json'])
        result.append({'marker':marker,'sample_id':sample,'fastq':str(dest),'reads':json.loads((folder/'qc.json').read_text())['reads']})
    atomic_json(out/'manifest.json',result);return result

def trim(cfg,merged):
    out=Path(cfg['_out'])/'trim';result=[]
    for row in merged:
        marker=row['marker'];sample=row['sample_id'];m=cfg['markers'][marker];src=Path(row['fastq']);folder=out/marker/sample;folder.mkdir(parents=True,exist_ok=True);dest=folder/(sample+'.fastq.gz')
        tool_version=subprocess.check_output(['cutadapt','--version'],text=True).strip() if m['trim']=='linked' else 'none'
        sig=digest({'algorithm':'trim-v1','code_sha256':sha(Path(__file__)),'input':sha(src),'parameters':m,'cutadapt_version':tool_version})
        if checked_checkpoint(folder,sig):print('RESUME trim',marker,sample,flush=True)
        else:
            temporary=folder/(sample+'.partial.fastq.gz')
            if m['trim']=='none':shutil.copyfile(src,temporary)
            else:
                F=m['forward'].upper();R=m['reverse'].upper();require(set(F+R)<=set('ACGTRYSWKMBDHVN'),'Primers must be verified IUPAC DNA')
                rc=lambda s:s.translate(str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN'))[::-1]
                require(min(len(F),len(R))>=m.get('overlap',16),'Primer shorter than requested overlap')
                cmd=['cutadapt','-j',str(cfg['threads']),'-e',str(m.get('error_rate',.15)),'-O',str(m.get('overlap',16)),'-g',f'F_R={F}...{rc(R)}','-g',f'R_F={R}...{rc(F)}','-m',str(m['min_length']),'-M',str(m['max_length']),'--untrimmed-output',str(folder/'untrimmed.fastq.gz'),'--too-short-output',str(folder/'too_short.fastq.gz'),'--too-long-output',str(folder/'too_long.fastq.gz'),'--json',str(folder/'cutadapt.json'),'-o',str(temporary),str(src)]
                command(cmd,folder/'cutadapt.log');atomic_json(folder/'command.json',cmd)
            n=bases=0;hist=collections.Counter()
            for h,s,p,q in fastq(temporary):n+=1;bases+=len(s);hist[len(s)]+=1
            temporary.replace(dest);atomic_json(folder/'qc.json',{'reads':n,'input_reads':row['reads'],'bases':bases,'technical_discarded':row['reads']-n,'primer_mode':m['trim'],'empty_sample':n==0});write_tsv(folder/'length_distribution.tsv',['length','reads'],[{'length':k,'reads':v} for k,v in sorted(hist.items())]);files=[f for f in folder.iterdir() if f.is_file() and f.name not in {'COMPLETE.json'} and '.partial' not in f.name];commit(folder,sig,files)
        qc=json.loads((folder/'qc.json').read_text());result.append(dict(marker=marker,sample_id=sample,fastq=str(dest),reads=qc['reads'],input_reads=row['reads']))
    atomic_json(out/'manifest.json',result);return result
