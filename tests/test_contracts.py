import contextlib,csv,gzip,io,json,os,sqlite3,tempfile,unittest
from pathlib import Path
import numpy as np
from pagleba.common import WorkflowError,fastq,atomic_json,write_tsv,read_tsv,sha,commit,checked_checkpoint,RANKS
from pagleba.classify import common_levels,classify_hits,mean_qscore
from pagleba.preprocess import rewrite_description,merge
from pagleba import config
from pagleba.functions import import_traits,run_functions
from pagleba.ncbi import prepare as ncbi_prepare,Taxonomy
from pagleba.statistics import prepare_jobs

class Contracts(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory(prefix='pagleba test ');self.root=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def cfg(self,profile='18s'):
        p=config.initialize(self.root/'project',profile);return config.load(p)
    def hit(self,ref,score=100,identity=.99,coverage=.99):return dict(query='q',ref=ref,score=score,identity=identity,coverage=coverage)
    def lineages(self):return {'a':['Eukaryota','Fungi','P','C','O','F','GenusA','SpeciesA'],'b':['Eukaryota','Fungi','P','C','O','F','GenusB','SpeciesB'],'c':['Eukaryota','Fungi','P','C','O','F','GenusA','SpeciesC']}
    def test_lca_same_genus_species_not_asserted(self):
        levels=common_levels([self.lineages()['a'],self.lineages()['c']]);self.assertEqual(levels[6],'GenusA');self.assertEqual(levels[7],'')
    def test_lca_disagreeing_genera_retains_family(self):self.assertEqual(common_levels([self.lineages()['a'],self.lineages()['b']])[5:8],['F','',''])
    def test_candidate_limit_withholds_genus(self):
        m=dict(identity=.95,coverage=.8,score_ratio=.99,secondary=1);r=classify_hits([self.hit('a'),self.hit('c')],self.lineages().__getitem__,m);self.assertTrue(r['CandidateLimitReached']);self.assertEqual(r['Status'],'eukaryotic_unresolved_genus')
    def test_invalid_high_score_does_not_override_valid_hit(self):
        m=dict(identity=.95,coverage=.8,score_ratio=.99,secondary=50);r=classify_hits([self.hit('a',300,identity=.8),self.hit('b',100)],self.lineages().__getitem__,m);self.assertEqual(r['BestReference'],'b')
    def test_no_acceptable_alignment_is_explicit(self):
        r=classify_hits([self.hit('a',coverage=.2)],self.lineages().__getitem__,dict(identity=.95,coverage=.8,score_ratio=.99));self.assertEqual(r['Status'],'no_acceptable_alignment')
    def test_plants_are_eukaryotic_assignments(self):
        r=classify_hits([self.hit('plant')],lambda x:['Eukaryota','Viridiplantae','Streptophyta','C','O','F','PlantGenus','S'],dict(identity=.95,coverage=.8,score_ratio=.99));self.assertEqual(r['Status'],'eukaryotic_genus');self.assertIn('Viridiplantae',r['LevelsJSON'])
    def test_uuid_sequence_and_quality_not_rewritten(self):
        h='@uuid_barcode01 barcode=barcode01 barcode_alias=barcode01';out,n=rewrite_description(h,'barcode01','barcode21');self.assertEqual(out,'@uuid_barcode01 barcode=barcode21 barcode_alias=barcode21');self.assertEqual(n,2)
    def test_no_broad_barcode_substring_substitution(self):
        h='@UUID something=notbarcode01x barcode=barcode01';out,n=rewrite_description(h,'barcode01','barcode09');self.assertIn('notbarcode01x',out);self.assertEqual(n,1)
    def test_fastq_duplicate_ids_preserved(self):
        p=self.root/'x.fastq.gz'
        with gzip.open(p,'wt') as f:f.write('@same\nACGT\n+\nIIII\n'*2)
        self.assertEqual(len(list(fastq(p))),2)
    def test_fastq_quality_mismatch_rejected(self):
        p=self.root/'bad.fastq';p.write_text('@x\nACGT\n+\nIII\n')
        with self.assertRaises(WorkflowError):list(fastq(p))
    def test_invalid_fastq_nucleotide_rejected(self):
        p=self.root/'invalid.fastq';p.write_text('@x\nAC Z\n+\nIIII\n')
        with self.assertRaises(WorkflowError):list(fastq(p))
    def test_its_does_not_inherit_18s_length_limits(self):
        c=self.cfg('its');self.assertIsNone(c['markers']['its']['min_length']);self.assertEqual(c['taxonomy_filter']['exclude_clades'],[])
    def test_cli_validation_error_is_reported_without_traceback(self):
        from pagleba.cli import main
        out=self.root/'existing.tsv';out.write_text('sample_id\nfixture\n')
        manifest=self.root/'samples.tsv';write_tsv(manifest,['sample_id'],[dict(sample_id='sample1')])
        err=io.StringIO()
        with contextlib.redirect_stderr(err):result=main(['metadata','--samples',str(manifest),'--out',str(out)])
        self.assertEqual(result,2);self.assertIn('Output exists',err.getvalue());self.assertNotIn('Traceback',err.getvalue())
    def test_generic_ssu_requires_explicit_assay_lengths(self):
        self.assertIsNone(config.PROFILES['18s']['min_length'])
        self.assertIsNone(config.PROFILES['18s']['max_length'])
        self.assertEqual(config.PROFILES['18s-f566-r1200']['min_length'],400)
        self.assertEqual(config.PROFILES['18s-f566-r1200']['max_length'],1200)
    def test_checkpoint_detects_corruption(self):
        d=self.root/'checkpoint';d.mkdir();f=d/'x';f.write_text('x');commit(d,'signature',[f]);self.assertTrue(checked_checkpoint(d,'signature'));f.write_text('corrupted')
        with self.assertRaises(WorkflowError):checked_checkpoint(d,'signature')
    def test_checkpoint_detects_scientific_parameter_changes(self):
        d=self.root/'checkpoint';d.mkdir();f=d/'x';f.write_text('x');commit(d,'old',[f])
        with self.assertRaises(WorkflowError):checked_checkpoint(d,'new')
    def test_reused_file_in_manifest_rejected(self):
        c=self.cfg();f=self.root/'r.fastq';f.write_text('@x barcode=barcode01\nACGT\n+\nIIII\n');rows=[dict(sample_id='s1',marker='ssu',run_id='r',barcode='barcode01',fastq=str(f)),dict(sample_id='s2',marker='ssu',run_id='r',barcode='barcode01',fastq=str(f))];write_tsv(config.resolve(c,c['samples']),list(rows[0]),rows)
        with self.assertRaises(WorkflowError):config.manifest(c)
    def test_merge_multiple_runs_counts_without_deduplication(self):
        c=self.cfg();rows=[]
        for i in range(2):
            f=self.root/f'{i}.fastq';f.write_text('@UUID barcode=barcode01\nACGT\n+\nIIII\n');rows.append(dict(sample_id='sampleA',marker='ssu',run_id=str(i),barcode='barcode01',header_barcode='barcode21',files=[f]))
        result=merge(c,rows);self.assertEqual(result[0]['reads'],2);records=list(fastq(result[0]['fastq']));self.assertTrue(all(h.startswith('@UUID ') for h,s,p,q in records));self.assertTrue(all('barcode21' in h for h,s,p,q in records));self.assertEqual(merge(c,rows),result)
    def test_qscore_averages_error_probabilities(self):self.assertAlmostEqual(mean_qscore('II'),40);self.assertLess(mean_qscore('!I'),4)
    def test_eukfunc_import_does_not_assign_unknown(self):
        f=self.root/'ref.tsv';write_tsv(f,['taxa','main_functional_class','level'],[dict(taxa='A',main_functional_class='unknown',level='genus'),dict(taxa='B',main_functional_class='predator',level='genus')]);out=self.root/'traits.tsv';import_traits('eukfunc',f,out);self.assertEqual(len(read_tsv(out)),1)
    def test_traits_do_not_multiply_repeated_reference_evidence(self):
        c=self.cfg();table=self.root/'counts.tsv';row=dict(Feature='PG_x',**dict(zip(RANKS,['Eukaryota','Fungi','P','C','O','F','G',''])),sample01=10);write_tsv(table,list(row),[row]);ann=self.root/'ann.tsv';a=dict(Source='test',Rank='Genus',Taxon='G',Trait='function',Weight=1,Evidence='fixture');write_tsv(ann,list(a),[a,a]);c['functions']=[dict(id='traits',marker='ssu',annotations=str(ann),method='long')];result=run_functions(c,[('ssu',table)]);counts=read_tsv(Path(result[0][1]).with_name('trait_counts.tsv'));self.assertEqual(float(counts[0]['sample01']),10)
    def test_unclassified_genus_not_used_for_function(self):
        c=self.cfg();table=self.root/'counts.tsv';row=dict(Feature='PG_x',**dict(zip(RANKS,['Eukaryota','Fungi','P','C','O','F','Unclassified_genus',''])),sample01=10);write_tsv(table,list(row),[row]);ann=self.root/'ann.tsv';a=dict(Source='test',Rank='Genus',Taxon='Unclassified_genus',Trait='function',Weight=1,Evidence='fixture');write_tsv(ann,list(a),[a]);c['functions']=[dict(id='traits',marker='ssu',annotations=str(ann),method='long')];run_functions(c,[('ssu',table)]);self.assertEqual(read_tsv(Path(c['_out'])/'functions/traits/trait_counts.tsv'),[])
    def test_singleton_group_distances_are_descriptive(self):
        c=self.cfg();write_tsv(config.resolve(c,c['metadata']),['sample_id','group'],[dict(sample_id='a',group='control'),dict(sample_id='b',group='treatment')]);atomic_json(config.resolve(c,c['tasks']),[dict(id='test',group='group',subset={},covariates=[],strata=None,blocks=['ssu'])]);table=self.root/'counts.tsv';base=dict(zip(RANKS,['Eukaryota','Fungi','P','C','O','F','G','']));write_tsv(table,['Feature']+RANKS+['a','b'],[dict(Feature='f1',**base,a=3,b=1),dict(Feature='f2',**base,a=1,b=3)]);c['statistics'].update(rarefaction_depth=2,rarefaction_repeats=2);jobs=prepare_jobs(c,[('ssu',table)],[]);self.assertEqual(jobs[0]['inferential_status'],'descriptive_only');summary=read_tsv(Path(jobs[0]['folder'])/'all_group_pairs.tsv');within=[r for r in summary if r['Group_A']==r['Group_B']];self.assertTrue(all(r['Mean']=='NA' for r in within));between=next(r for r in summary if r['Group_A']!=r['Group_B']);self.assertAlmostEqual(float(between['Mean']),.5)
    def test_ncbi_synonym_and_genus_split_are_not_conflated(self):
        nodes=self.root/'nodes.dmp';names=self.root/'names.dmp';nodes.write_text('1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tdomain\t|\n10\t|\t2\t|\tgenus\t|\n11\t|\t2\t|\tgenus\t|\n');names.write_text('1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\t\t|\tscientific name\t|\n10\t|\tCurrentgenus\t|\t\t|\tscientific name\t|\n10\t|\tOldsynonym\t|\t\t|\tsynonym\t|\n11\t|\tSplitgenus\t|\t\t|\tscientific name\t|\n');db=self.root/'ncbi.sqlite';ncbi_prepare(nodes,names,db);tx=Taxonomy(db);self.assertEqual(tx.resolve('Oldsynonym','Genus')[0][0],10);self.assertEqual(tx.resolve('Splitgenus','Genus')[0][0],11)

if __name__=='__main__':unittest.main()
