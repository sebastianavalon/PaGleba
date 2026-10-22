import json,os,shutil,tempfile,unittest
from pathlib import Path
from pagleba.common import read_tsv,WorkflowError,sha
from pagleba.demo import create
from pagleba.config import load
from pagleba.workflow import run
from pagleba.classify import classify

@unittest.skipUnless(os.environ.get('PAGLEBA_INTEGRATION')=='1' and shutil.which('minimap2') and shutil.which('Rscript'),'Set PAGLEBA_INTEGRATION=1 in the full supported environment')
class Integration(unittest.TestCase):
    def test_two_markers_pilot_resume_models_plants_and_small_groups(self):
        with tempfile.TemporaryDirectory(prefix='PaGleba integration ') as tmp:
            cfgpath=create(Path(tmp)/'project',12,2);cfg=load(cfgpath);cfg['statistics']['network_bootstrap']=40
            annotation=cfgpath.parent/'empty_annotations.tsv';annotation.write_text('Source\tRank\tTaxon\tTrait\tWeight\tEvidence\nfixture\tGenus\tNoSuchGenus\tno_match\t1\ttest\n');cfg['functions'].append(dict(id='empty_annotations',marker='ssu',annotations=str(annotation),method='long'));tasks=json.loads((cfgpath.parent/'tasks.json').read_text());tasks[0]['blocks'].append('empty_annotations');(cfgpath.parent/'tasks.json').write_text(json.dumps(tasks))
            run(cfg,pilot=1);out=Path(cfg['_out']);self.assertTrue((out/'taxonomy/pilot.json').exists());self.assertFalse((out/'taxonomy/CLASSIFICATION_COMPLETE.json').exists())
            run(cfg,trees=True)
            for marker in ['ssu','its']:
                table=read_tsv(out/'tables'/marker/'otu_genus_counts.tsv');self.assertTrue(any(r['Kingdom']=='Viridiplantae' for r in table));samples=[f'sample{i:02d}' for i in range(1,13)];self.assertEqual(sum(int(r[s]) for r in table for s in samples),12*24);self.assertTrue(json.loads((out/'tables'/marker/'count_reconciliation.json').read_text())['passed'])
            stats=read_tsv(out/'statistics/PERMANOVA_summary.tsv');self.assertTrue(all(r['Status']=='descriptive_only' for r in stats if r['Task']=='small_subset'));self.assertTrue(all(r['P']=='NA' for r in stats if r['Task']=='small_subset'))
            self.assertEqual(next(r for r in stats if r['Block']=='empty_annotations')['Status'],'not_estimable_empty_block')
            models=read_tsv(out/'models/MODEL_STATUS.tsv');self.assertIn('not_inferred',next(r for r in models if r['Task']=='small_subset')['Network']);self.assertTrue((out/'models/all_groups/PLS_raw/paths.tsv').exists());self.assertTrue((out/'models/cross_block_concordance.tsv').exists());self.assertTrue((out/'report/REPORT.html').exists());self.assertIn('SYNTHETIC TEST DATA',(out/'report/REPORT.html').read_text());self.assertEqual(json.loads((out/'report/validation.json').read_text())['missing_assets'],[])
            before={str(p.relative_to(out)):sha(p) for p in (out/'tables').rglob('*.tsv')};run(cfg,stages=['merge','trim','classify','tables','functions','statistics','models']);after={str(p.relative_to(out)):sha(p) for p in (out/'tables').rglob('*.tsv')};self.assertEqual(before,after)
            changed=load(cfgpath);changed['markers']['its']['database']=changed['markers']['ssu']['database']
            with self.assertRaises(WorkflowError):classify(changed,json.loads((out/'trim/manifest.json').read_text()))
