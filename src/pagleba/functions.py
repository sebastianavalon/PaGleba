"""Audited positive taxon-to-trait associations; missing annotations are unknown."""
import collections,json,re,sqlite3
from pathlib import Path
import numpy as np
from .common import RANKS,read_tsv,write_tsv,require,sha,digest,atomic_json,checked_checkpoint,commit
from .config import resolve

ANNOTATION_FIELDS=['Source','Rank','Taxon','Trait','Weight','Evidence','ContextFamily','ContextPhylum','ContextKingdom']
def clean(s):return ' '.join((s or '').replace('_',' ').split()).casefold()
def usable(s):return bool(s) and not any(t in clean(s) for t in ['unclassified','incertae','uncultured','unknown','unidentified']) and not bool(re.fullmatch(r'[\d_.-]+',s))
def import_traits(fmt,reference,out,fields=None):
    rows=read_tsv(reference);result=[]
    if fmt=='eukfunc':
        for r in rows:
            require({'taxa','main_functional_class','level'}<=set(r),'Expected EukFunc DBc_minimal/DBc_main/DBc TSV export')
            if r['main_functional_class'] in {'','unknown','unassigned'}:continue
            rank=r['level'].title();rank=rank if rank in RANKS else 'Any'
            result.append(dict(Source='EukFunc',Rank=rank,Taxon=r['taxa'],Trait='main_functional_class::'+r['main_functional_class'],Weight=1,Evidence='Database clade annotation; not measured function',ContextFamily=''))
    elif fmt=='fungaltraits':
        for r in rows:
            genus=r.get('GENUS') or r.get('Genus');require(genus is not None,'Expected FungalTraits genus TSV (GENUS column)')
            for field in fields or ['primary_lifestyle','Secondary_lifestyle','Endophytic_interaction_capability_template','Plant_pathogenic_capacity_template','Decay_substrate_template','Decay_type_template','Aquatic_habitat_template','Growth_form_template']:
                for val in re.split(r',\s*(?![^()]*\))',r.get(field,'')):
                    if val.strip().casefold() in {'','unknown','unspecified','na','n/a','none_available'}:continue
                    result.append(dict(Source='FungalTraits',Rank='Genus',Taxon=genus,Trait=field+'::'+val.strip(),Weight=1,Evidence='Reference genus observation; unknown is not a negative phenotype',ContextFamily=r.get('Family',''),ContextPhylum=r.get('Phylum',''),ContextKingdom=''))
    else:
        require(rows and {'Source','Rank','Taxon','Trait'}<=set(rows[0]),'Long TSV requires Source,Rank,Taxon,Trait')
        for r in rows:result.append(dict(r,Weight=r.get('Weight') or 1,Evidence=r.get('Evidence','User-supplied positive association'),ContextFamily=r.get('ContextFamily','')))
    write_tsv(out,ANNOTATION_FIELDS,result);atomic_json(str(out)+'.source.json',{'format':fmt,'reference_sha256':sha(reference),'rows':len(result),'method':'Explicit exact taxon/rank annotations; no fuzzy name matching'})

def run_functions(cfg,tables):
    available=dict(tables);outputs=[]
    for spec in cfg.get('functions',[]):
        name=spec['id'];require(re.fullmatch(r'[A-Za-z0-9_.-]+',name),'Unsafe functional block ID');marker=spec['marker'];require(marker in available,f'No table for function marker {marker}');source=available[marker];annotation=resolve(cfg,spec['annotations']);folder=Path(cfg['_out'])/'functions'/name;folder.mkdir(parents=True,exist_ok=True)
        aliasfile=resolve(cfg,spec['aliases']) if spec.get('aliases') else None
        sig=digest({'algorithm':'traits-v1','code_sha256':sha(Path(__file__)),'counts':sha(source),'annotations':sha(annotation),'aliases':sha(aliasfile) if aliasfile else None,'spec':spec})
        if checked_checkpoint(folder,sig):outputs.append((name,folder/'trait_fraction.tsv'));continue
        rows=read_tsv(source);require(rows,'Empty taxonomy table for functional analysis');samples=[k for k in rows[0] if k not in ['Feature']+RANKS];totals=np.sum([[float(r[s]) for s in samples] for r in rows],axis=0)
        annotations=read_tsv(annotation);index=collections.defaultdict(list)
        for a in annotations:
            require({'Source','Rank','Taxon','Trait'}<=set(a),'Invalid annotation schema');rank=a['Rank'].title();require(rank in RANKS+['Any','Feature'],'Unknown annotation rank');require(0<float(a.get('Weight') or 1)<=1,'Positive trait weights must be in (0,1]');index[(rank,clean(a['Taxon']))].append(a)
        aliases={};audit=[]
        if aliasfile:
            for a in read_tsv(aliasfile):
                require({'Rank','Original','Canonical','Evidence'}<=set(a) and a['Evidence'],'Alias file needs Rank,Original,Canonical,Evidence; only verified same-taxID aliases')
                key=(a['Rank'].title(),clean(a['Original']));require(key not in aliases,'Conflicting alias records');aliases[key]=a
        if aliases:
            rebuilt=collections.defaultdict(list)
            for (rank,taxon),annrows in index.items():
                alias=aliases.get((rank,taxon));rebuilt[(rank,clean(alias['Canonical']) if alias else taxon)].extend(annrows)
            index=rebuilt
        membership=collections.defaultdict(dict);feature_rows=[];matched=set()
        for r in rows:
            candidates=[(8,ann,r['Feature'],r['Feature']) for ann in index.get(('Feature',clean(r['Feature'])),[])]
            for rank in RANKS[:-1]:
                if not usable(r[rank]):continue
                raw=r[rank];a=aliases.get((rank,clean(raw)));taxon=a['Canonical'] if a else raw
                if a:audit.append(dict(Feature=r['Feature'],Rank=rank,Original=raw,Canonical=taxon,Evidence=a['Evidence']))
                for ann in index.get((rank,clean(taxon)),[])+index.get(('Any',clean(taxon)),[]):
                    if ann['Source']=='FungalTraits':
                        if r.get('Domain')!='Eukaryota':continue
                        phylum=ann.get('ContextPhylum','');input_phylum=r.get('Phylum','')
                        if phylum and usable(input_phylum) and clean(phylum)!=clean(input_phylum):continue
                        if clean(r.get('Kingdom',''))!='fungi' and not(phylum and usable(input_phylum) and clean(phylum)==clean(input_phylum)):continue
                    if ann['Source']=='EukFunc' and r.get('Domain')!='Eukaryota':continue
                    if ann['Source'].startswith('MACADAM') and r.get('Domain') not in {'Bacteria','Archaea'}:continue
                    family=ann.get('ContextFamily','')
                    if family and usable(r.get('Family','')) and clean(family)!=clean(r['Family']):continue
                    candidates.append((RANKS.index(rank),ann,raw,taxon))
            # EukFunc is a most-specific-clade classifier; do not accumulate contradictory parent and child primary classes.
            if spec.get('method')=='eukfunc' and candidates:
                deepest=max(v[0] for v in candidates);candidates=[v for v in candidates if v[0]==deepest]
                classes={v[1]['Trait'] for v in candidates};require(len(classes)<=1,f'Ambiguous EukFunc clade for {r["Feature"]}; inspect reference aliases')
            for _,a,raw,taxon in candidates:
                key=(a['Source'],a['Trait']);weight=float(a.get('Weight') or 1);membership[key][r['Feature']]=max(weight,membership[key].get(r['Feature'],0));matched.add(r['Feature']);feature_rows.append(dict(Feature=r['Feature'],Source=a['Source'],Trait=a['Trait'],Input_Taxon=raw,Matched_Taxon=taxon,Rank=a['Rank'],Weight=weight,Evidence=a.get('Evidence','')))
        byid={r['Feature']:r for r in rows};rawrows=[];fraction=[]
        for (src,trait),features in sorted(membership.items()):
            vector=sum((np.array([float(byid[f][s]) for s in samples])*w for f,w in features.items()),np.zeros(len(samples)));head=dict(Feature='FN_'+digest((src,trait))[:16],Source=src,Trait=trait);rawrows.append(dict(head,**dict(zip(samples,vector))));fraction.append(dict(head,**{s:(float(vector[j]/totals[j]) if totals[j]>0 else 'NA') for j,s in enumerate(samples)}))
        fields=['Feature','Source','Trait']+samples;write_tsv(folder/'trait_counts.tsv',fields,rawrows);write_tsv(folder/'trait_fraction.tsv',fields,fraction);write_tsv(folder/'trait_percent.tsv',fields,[dict(r,**{s:100*r[s] if r[s]!='NA' else 'NA' for s in samples}) for r in fraction]);write_tsv(folder/'feature_trait_evidence.tsv',['Feature','Source','Trait','Input_Taxon','Matched_Taxon','Rank','Weight','Evidence'],feature_rows);write_tsv(folder/'name_alias_audit.tsv',['Feature','Rank','Original','Canonical','Evidence'],audit)
        covered=sum((np.array([float(byid[f][s]) for s in samples]) for f in matched),np.zeros(len(samples)));write_tsv(folder/'annotation_coverage.tsv',['Sample','Total_Target_Reads','Reads_With_Any_Annotation','Coverage_Fraction'],[dict(Sample=s,Total_Target_Reads=totals[j],Reads_With_Any_Annotation=covered[j],Coverage_Fraction=covered[j]/totals[j] if totals[j] else 'NA') for j,s in enumerate(samples)])
        atomic_json(folder/'interpretation.json',{'unmatched_taxa':'unknown, not functionally negative','overlapping_traits':True,'denominator':'all retained target-domain reads before functional assignment','taxonomic_filter':'only explicit project taxonomy filter; no automatic plant removal','inference':'reference annotation potential; not measured genes, expression, activity or causation','matched_rows':len(matched),'input_rows':len(rows)})
        commit(folder,sig,[f for f in folder.iterdir() if f.is_file() and f.name!='COMPLETE.json']);outputs.append((name,folder/'trait_fraction.tsv'))
    return outputs
