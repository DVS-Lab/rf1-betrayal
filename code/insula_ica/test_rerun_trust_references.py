"""Verify stock-query scoping, anatomical isolation, and costly-job guards."""
from contextlib import ExitStack
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import rerun_trust_references as r


class StockRerunTests(unittest.TestCase):
    def fixture(self, root):
        raw=root/'project/bids';session=raw/'sub-10001/ses-01'
        for name in ('anat','func','fmap'):(session/name).mkdir(parents=True)
        (raw/'dataset_description.json').write_text('{}')
        (session/'anat/sub-10001_ses-01_T1w.nii.gz').write_bytes(b'raw anatomy')
        files={}
        for run in (1,2):
            for echo in (1,2,3,4):
                stem=f'sub-10001_ses-01_task-trust_run-{run}_echo-{echo}_part-mag'
                for suffix in ('_bold.nii.gz','_bold.json','_sbref.nii.gz','_sbref.json'):
                    path=session/'func'/(stem+suffix);path.write_bytes(b'raw fixture')
                    if run==2 and suffix=='_bold.nii.gz':files[str(path)]={}
        old=root/'old/fmriprep-24/sub-10001/ses-01'
        (old/'anat').mkdir(parents=True)
        for name in ('sub-10001_ses-01_desc-preproc_T1w.nii.gz','sub-10001_ses-01_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5'):
            (old/'anat'/name).write_bytes(b'cached anatomy')
        stock=raw.parent/'code/fmriprep_config.json';stock.parent.mkdir()
        stock.write_bytes(b'{\n "sbref":{"datatype":"func","suffix":"sbref","part":[null,"mag"]},\n "bold":{"datatype":"func","suffix":"bold","part":[null,"mag"]}\n}\n')
        (stock.parent/'fmriprep.sh').write_text('# original wrapper, unchanged\n')
        conf=root/'original-confounds.tsv';conf.write_text('1\n2\n')
        row=dict(key='sub-10001_ses-01_task-trust_run-2',subject='sub-10001',run=2,raw_images=files,
                 new_bold=str(old/'func/sub-10001_ses-01_task-trust_run-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz'),
                 confounds=str(conf),confounds_sha256=r.core.sha(conf))
        return raw,stock,row

    def test_view_preserves_stock_filter_bytes_and_natural_sbref_matching(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);raw,stock,row=self.fixture(root)
            before={str(p):p.read_bytes() for p in raw.rglob('*') if p.is_file()}
            case=root/'case';anat,identities=r.stage_inputs(row,case,raw,stock)
            self.assertEqual((case/'fmriprep_config.json').read_bytes(),stock.read_bytes())
            func=case/'bids-view/sub-10001/ses-01/func'
            self.assertEqual(len(list(func.glob('*_bold.nii.gz'))),4)
            self.assertTrue(all('_run-2_' in p.name for p in func.glob('*_bold.nii.gz')))
            self.assertEqual(len(list(func.glob('*_sbref.nii.gz'))),8)
            self.assertNotIn('run',json.loads((case/'fmriprep_config.json').read_text())['sbref'])
            cache=case/'anatomy-cache'
            self.assertFalse(list(cache.rglob('func')))
            self.assertTrue(list(cache.rglob('*.h5')))
            self.assertTrue(all(Path(x['path']).is_relative_to(anat) for x in identities))
            self.assertEqual(before,{str(p):p.read_bytes() for p in raw.rglob('*') if p.is_file()})
            with self.assertRaises(FileExistsError):r.stage_inputs(row,case,raw,stock)

    def test_command_reuses_only_anatomy_and_does_not_override_reference_queries(self):
        argv=r.command('apptainer',Path('/fresh'),Path('/work'),'sub-10001',Path('/old/anat'))
        self.assertEqual(argv[argv.index('--derivatives')+1],'anatomy=/rebuild/anatomy-cache')
        self.assertNotIn('--ignore',argv)
        self.assertNotIn('--random-seed',argv)
        self.assertEqual(argv[argv.index('--bids-filter-file')+1],'/rebuild/fmriprep_config.json')
        self.assertIn('/old/anat:/old/anat:ro',argv)
        self.assertIn('fmriprep-24.1.1.simg',str(r.IMAGE))
        compile(r.ANATOMY_PREFLIGHT,'installed-anatomy-preflight','exec')

    def test_repair_selection_covers_all_108_remaining_and_rejects_unknowns(self):
        rows=[dict(key=f'sub-{10000+i}_ses-01_task-trust_run-2',subject=f'sub-{10000+i}',run=2) for i in range(108)]
        rows.append(dict(key='sub-10834_ses-01_task-trust_run-2',subject='sub-10834',run=2))
        audit=[dict(key=x['key'],status='different_run_reference') for x in rows]
        repaired={'sub-10834_ses-01_task-trust_run-2'}
        self.assertEqual(len(r.choose_runs(rows,audit,repaired)),108)
        audit[0]['status']='unknown_no_retained_log_evidence'
        with self.assertRaisesRegex(ValueError,'unknown'):r.choose_runs(rows,audit,repaired)

    def test_actual_child_is_stopped_if_anatomical_normalization_is_requested(self):
        with tempfile.TemporaryDirectory() as directory:
            log=Path(directory)/'run.log'
            with self.assertRaisesRegex(ValueError,'new anatomical normalization'):
                r.execute([sys.executable,'-u','-c',
                    "import time; print('ANAT Stage 4: Preparing normalization workflow'); time.sleep(5)"],dict(os.environ),log)
            self.assertIn('Preparing normalization',log.read_text())
            r.execute([sys.executable,'-c',"print('ANAT Stage 4: Found pre-computed registrations')"],dict(os.environ),log)

    def test_absent_sbref_uses_the_stock_bold_derived_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);key='sub-10001_ses-01_task-trust_run-2'
            row=dict(subject='sub-10001',key=key,new_bold=str(root/'fmriprep-24/sub-10001/ses-01/func/bold.nii.gz'))
            (root/'sub-10001.fmriprep.log').write_text(f'No single-band-reference found for {key}_echo-1_part-mag_bold.nii.gz.\n')
            self.assertEqual(r.sbref_selection([row])[0]['status'],'bold_derived_reference_no_sbref')

    def test_successful_pipeline_records_actual_reference_and_keeps_anatomical_inputs(self):
        with tempfile.TemporaryDirectory() as directory,ExitStack() as p:
            root=Path(directory);raw,stock,row=self.fixture(root)
            source=root/'code/insula_ica/records'/r.SOURCE;source.mkdir(parents=True)
            (source/'input-registry.json').write_text(json.dumps([row]))
            image=root/'fmriprep-24.1.1.simg';image.write_bytes(b'image fixture')
            (source/'container.json').write_text(json.dumps(dict(sha256=r.rebuild.sha(image))))
            repair=root/'code/insula_ica/records/linux2-registration-fix-sub-10834-20261028T151154Z-acdf91';repair.mkdir()
            (repair/'review-decision.json').write_text(json.dumps(dict(assessment='accept_repaired_pair_for_sub_10834_analysis_use',selected_inputs=[])))
            (repair/'cohort-sbref-selection-audit.json').write_text(json.dumps([dict(key=row['key'],status='different_run_reference')]))
            record=root/'record';record.mkdir()
            p.enter_context(patch.object(r.core,'ROOT',root));p.enter_context(patch.object(r,'RAW',raw))
            p.enter_context(patch.object(r,'WORK',root/'work'));p.enter_context(patch.object(r,'STOCK',stock))
            p.enter_context(patch.object(r,'IMAGE',image));p.enter_context(patch.dict(os.environ,FSLDIR='/fixture/fsl'))
            p.enter_context(patch.object(r,'choose_runs',return_value=[row]))
            p.enter_context(patch.object(r.shutil,'which',side_effect=lambda x:'/usr/bin/'+x))
            case=root/'derivatives/insula_ica/trust-stock-filter-fixture'/row['subject']
            def preflight(*args,**kwargs):
                (case/'anatomy-preflight.json').write_text(json.dumps(dict(anatomical_registration_nodes_absent=True)))
            p.enter_context(patch.object(r.subprocess,'run',side_effect=preflight))
            def execute(argv,env,log):
                log.write_text('ANAT Stage 4: Found pre-computed registrations\n'+
                    f"Report - setting before (/rebuild/fmriprep-24/{row['subject']}/ses-01/func/{row['key']}_part-mag_desc-coreg_boldref.nii.gz) "
                    f"and after (/bids/{row['subject']}/ses-01/func/{row['key']}_echo-1_part-mag_sbref.nii.gz) images\n")
            p.enter_context(patch.object(r,'execute',side_effect=execute))
            p.enter_context(patch.object(r.rebuild,'verify_outputs'))
            p.enter_context(patch.object(r.rebuild,'warpkit_inputs',return_value=dict(ready=True)))
            fit=p.enter_context(patch.object(r.rebuild,'fit_run'))
            p.enter_context(patch.object(r,'qc_run',return_value=dict(review_flags=[])))
            p.enter_context(patch.object(r,'montage',side_effect=lambda q,path,fsldir:path.write_bytes(b'png fixture')))
            state=dict(subjects={});r.run(record,state,'fixture')
            self.assertEqual(state['state'],'functional_reruns_complete_qc_pending')
            self.assertEqual(state['completed'],1);self.assertEqual(state['failed'],0)
            result=state['subjects'][row['subject']]
            self.assertEqual(result['reference_selection']['status'],'same_run')
            self.assertTrue(result['anatomical_reuse_verified'])
            self.assertEqual(fit.call_args.kwargs,dict(strict_historical=False))
            self.assertEqual((record/'fmriprep_config.json').read_bytes(),stock.read_bytes())
            self.assertEqual((Path(row['new_bold']).parent.parent/'anat/sub-10001_ses-01_desc-preproc_T1w.nii.gz').read_bytes(),b'cached anatomy')


if __name__=='__main__':unittest.main()
