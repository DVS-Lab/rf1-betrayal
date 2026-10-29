from contextlib import ExitStack, nullcontext, redirect_stdout
import io
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

import recover_ugr_unsmoothed as r
import test_unsmoothed_feat as feat_fixture
import test_rerun_trust_references as stock_fixture


class UGRRecoveryTests(unittest.TestCase):
    def test_command_retains_stock_and_cached_anatomy_and_changes_only_task_resources(self):
        original=r.stock.command('apptainer',Path('/case'),Path('/scratch'),'sub-1',Path('/anatomy'))
        new=r.build_command('apptainer',Path('/case'),Path('/scratch'),'sub-1',Path('/anatomy'))
        for flag,value in (('--task-id','ugr'),('--nthreads','8'),('--omp-nthreads','4'),('--mem-mb','12000')):
            self.assertEqual(new[new.index(flag)+1],value)
            original[original.index(flag)+1]=value
        self.assertEqual(new,original)
        self.assertIn('anatomy=/rebuild/anatomy-cache',new)

    def test_UGR_view_keeps_only_UGR_sbrefs_with_unchanged_filter(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);raw,stock,row=stock_fixture.StockRerunTests().fixture(root)
            func=raw/'sub-10001/ses-01/func'
            for path in list(func.iterdir()):
                (func/path.name.replace('task-trust','task-ugr')).write_bytes(path.read_bytes())
            row['raw_images']={str(func/Path(p).name.replace('task-trust','task-ugr')):{} for p in row['raw_images']}
            case=root/'case'
            r.stock.stage_inputs(row,case,raw,stock,task='ugr')
            view=case/'bids-view/sub-10001/ses-01/func'
            self.assertEqual(len(list(view.glob('*_bold.nii.gz'))),4)
            self.assertEqual(len(list(view.glob('*_sbref.nii.gz'))),8)
            self.assertTrue(all('_task-ugr_' in p.name for p in view.glob('*_sbref.nii.gz')))
            self.assertEqual((case/'fmriprep_config.json').read_bytes(),stock.read_bytes())

    def test_UGR_Warpkit_checks_task_specific_BIDS_maps_and_associations(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);raw=root/'bids';sub='sub-1';key=sub+'_ses-01_task-ugr_run-1'
            func=raw/sub/'ses-01/func';func.mkdir(parents=True)
            fmap=raw/sub/'ses-01/fmap';fmap.mkdir()
            for echo in range(1,5):
                stem=key+f'_echo-{echo}_part-mag_bold'
                nib.save(nib.Nifti1Image(np.zeros((2,2,2,6),dtype='float32'),np.eye(4)),func/(stem+'.nii.gz'))
                (func/(stem+'.json')).write_text(json.dumps(dict(RepetitionTime=1.615,B0FieldSource='ugr1')))
            stem=sub+'_ses-01_acq-ugr_run-1'
            for suffix in ('fieldmap','magnitude'):
                nib.save(nib.Nifti1Image(np.ones((2,2,2),dtype='float32'),np.eye(4)),fmap/(stem+'_'+suffix+'.nii.gz'))
                (fmap/(stem+'_'+suffix+'.json')).write_text(json.dumps(dict(Units='Hz',B0FieldIdentifier='ugr1')))
            row=dict(key=key,subject='1',run=1,nvol=6,tr=1.615)
            audit=r.raw_inputs(row,raw)
            self.assertTrue(audit['warpkit']['ready']);self.assertEqual(len(audit['images']),4)
            self.assertIn('task-ugr',audit['warpkit']['retained_4d_fieldmaps'])

    def test_explicit_override_preserves_original_and_missing_input_FSFs(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row,text,_=feat_fixture.UnsmoothedTests().fixture(root)
            row['task']='ugr'
            with patch.object(r.sm0.g.residuals,'validate_model'):
                missing=r.sm0.plan(row,root)
                old_fsf=Path(missing['unsmoothed_fsf']);frozen=old_fsf.read_bytes()
                bold=root/'regenerated.nii.gz'
                nib.save(nib.Nifti1Image(np.ones((5,5,5,6),dtype='float32'),np.eye(4)),bold)
                override=dict(path=str(bold),identity=r.core.identity(bold,False),source_record='/explicit/recovery')
                ready=r.sm0.plan(row,root,override)
                self.assertEqual(ready['status'],'ready')
                self.assertNotEqual(ready['unsmoothed_fsf'],str(old_fsf))
                self.assertEqual(old_fsf.read_bytes(),frozen)
                self.assertEqual(Path(row['source_fsf']).read_text(),text)
                override['identity']['size']+=1
                with self.assertRaisesRegex(ValueError,'BOLD changed'):r.sm0.plan(row,root,override)

    def test_BOLD_QC_checks_actual_regenerated_image_mask_and_saved_confounds(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);old=root/'old.feat';old.mkdir();case=root/'case'
            key='sub-1_ses-01_task-ugr_run-1'
            name=key+'_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz'
            folder=case/'fmriprep-24/sub-1/ses-01/func';folder.mkdir(parents=True)
            rng=np.random.default_rng(3);base=rng.uniform(100,200,size=(5,5,5)).astype('float32')
            data=base[...,None]+rng.normal(size=(5,5,5,6)).astype('float32')
            image=nib.Nifti1Image(data,np.eye(4));image.header.set_zooms((1,1,1,1.615));nib.save(image,folder/name)
            for path,values in [(old/'mean_func.nii.gz',base),(old/'example_func.nii.gz',base),
                                (old/'mask.nii.gz',np.ones(base.shape,dtype='uint8')),
                                (folder/name.replace('_desc-preproc_bold','_desc-brain_mask'),np.ones(base.shape,dtype='uint8'))]:
                nib.save(nib.Nifti1Image(values,np.eye(4)),path)
            (old/'design.fsf').write_text(f'set feat_files(1) "/original/{name}"\n')
            conf=root/'confounds.txt';conf.write_text('1\n2\n')
            row=dict(subject='1',key=key,nvol=6,tr=1.615,active_feat=str(old),confounds=str(conf),confounds_sha256=r.core.sha(conf))
            masks={roi:np.ones(base.shape,bool) for roi in ('left','right','bilateral')}
            with patch.dict(os.environ,FSLDIR='/fixture'),patch.object(r.core,'anatomical_masks',return_value=(masks,{},'insula')):
                qc=r.bold_qc(row,case)
                self.assertEqual(qc['coverage'],dict(left=1.,right=1.,bilateral=1.))
                self.assertEqual(qc['review_flags'],[]);self.assertEqual(qc['identity'],r.core.identity(folder/name,False))
                data[0,0,0,0]=np.nan;nib.save(nib.Nifti1Image(data,image.affine,image.header),folder/name)
                with self.assertRaisesRegex(ValueError,'Nonfinite'):r.bold_qc(row,case)

    def workflow(self, directory, fail_pilot=False):
        root=Path(directory);records=root/'code/insula_ica/records';records.mkdir(parents=True)
        source=records/'prepared';source.mkdir()
        rows=[]
        for index in range(132):
            subject=str(10000+index)
            for run in range(1,2 if index<15 else 3):
                rows.append(dict(subject=subject,task='ugr',run=run,key=f'sub-{subject}_ses-01_task-ugr_run-{run}'))
        self.assertEqual(len(rows),249)
        (source/'accepted-inputs.json').write_text(json.dumps(rows))
        old=records/'old-sm0';old.mkdir();featroot=root/'derivatives/insula_ica/sm-0-feat-old';featroot.mkdir(parents=True)
        sentinel=featroot/'completed-Trust.txt';sentinel.write_bytes(b'previously completed Trust output')
        (old/'state.json').write_text(json.dumps(dict(output_root=str(featroot),source_record=str(source))))
        trust=records/'linux2-trust-regeneration-20261027T162428Z-ef53c5';trust.mkdir()
        (trust/'input-registry.json').write_text(json.dumps([dict(subject='sub-'+s,new_bold='/old/anatomy/bold') for s in sorted({r['subject'] for r in rows})]))
        image=root/'container.simg';image.write_bytes(b'container')
        (trust/'container.json').write_text(json.dumps(dict(sha256=r.stock.rebuild.sha(image))))
        raw=root/'upstream/bids';raw.mkdir(parents=True);(raw.parent/'code').mkdir()
        stock=raw.parent/'code/fmriprep_config.json';stock.write_text('{}')
        (stock.parent/'fmriprep.sh').write_text('# unchanged wrapper\n')
        output=root/'derivatives/insula_ica/ugr-anatomy-reuse-test';output.mkdir()
        record=records/'linux2-ugr-anatomy-reuse-test';record.mkdir()
        args=SimpleNamespace(sm0_record=old,fmriprep_jobs=32,feat_jobs=64,publish=True)
        state=dict(subjects={},runs={},feat_runs={},code_commit='test')
        def stage(row,case,**kwargs):
            self.assertEqual(kwargs['task'],'ugr');self.assertEqual(len(row['raw_images']),4)
            case.mkdir(parents=True);(case/'fmriprep_config.json').write_text('{}')
            return Path('/old/anatomy'),[]
        def preflight(argv,**kwargs):
            case=Path(next(a[:-len(':/rebuild')] for a in argv if a.endswith(':/rebuild')))
            (case/'anatomy-preflight.json').write_text(json.dumps(dict(anatomical_registration_nodes_absent=True)))
        def execute(argv,env,log):
            log.write_text('ANAT Stage 4: Found pre-computed registrations\nfMRIPrep finished successfully!\n')
        def qc(row,case):
            return dict(key=row['key'],path=str(case/'bold.nii.gz'),identity={},coverage=dict(left=.5 if fail_pilot else 1.,right=1.,bilateral=1.),
                        historical_frame=dict(spatial_r=.99),review_flags=[])
        def fit(item,out,rec):
            (rec/'qc'/(item['key']+'.json')).write_text(json.dumps(dict(coverage=dict(left=1.,right=1.,bilateral=1.))))
            return dict(state='complete')
        def combined(args,out,rec,s):
            self.assertEqual(out,featroot)
            s.update(state='sm0_FEAT_complete_QC_pending',runs={str(i):dict(state='complete') for i in range(490)})
            r.save(rec/'state.json',s)
        with ExitStack() as p,redirect_stdout(io.StringIO()):
            for name,value in [('ROOT',root),('RECORDS',records),('WORK',root/'scratch')]:p.enter_context(patch.object(r,name,value))
            p.enter_context(patch.object(r.stock,'RAW',raw));p.enter_context(patch.object(r.stock,'STOCK',stock));p.enter_context(patch.object(r.stock,'IMAGE',image))
            p.enter_context(patch.dict(os.environ,FSLDIR='/fixture/fsl'))
            p.enter_context(patch.object(r.shutil,'which',side_effect=lambda name:'/bin/'+name))
            p.enter_context(patch.object(r,'available_gib',return_value=128))
            p.enter_context(patch.object(r.Budget,'reserve',side_effect=lambda *a:nullcontext()))
            p.enter_context(patch.object(r,'raw_inputs',side_effect=lambda row,raw:dict(images=['echo']*4,warpkit=dict(ready=True))))
            p.enter_context(patch.object(r.stock,'stage_inputs',side_effect=stage))
            p.enter_context(patch.object(r.subprocess,'run',side_effect=preflight))
            execute_mock=p.enter_context(patch.object(r.stock,'execute',side_effect=execute))
            p.enter_context(patch.object(r,'bold_qc',side_effect=qc))
            p.enter_context(patch.object(r,'sbref_selection',return_value=[dict(status='same_run')]))
            p.enter_context(patch.object(r.sm0,'plan',side_effect=lambda row,*a:dict(row,status='ready')))
            feat_mock=p.enter_context(patch.object(r.sm0,'fit',side_effect=fit))
            combined_mock=p.enter_context(patch.object(r.sm0,'run',side_effect=combined))
            publishing=p.enter_context(patch.object(r,'publish'))
            if fail_pilot:
                with self.assertRaisesRegex(ValueError,'pilot coverage/alignment'):r.run(args,record,state,output)
                self.assertEqual(execute_mock.call_count,3)
                feat_mock.assert_not_called();combined_mock.assert_not_called();publishing.assert_not_called()
            else:
                r.run(args,record,state,output)
                self.assertEqual(execute_mock.call_count,249)
                self.assertEqual(feat_mock.call_count,249)
                self.assertEqual(len(state['subjects']),132)
                self.assertEqual(state['state'],'UGR_and_sm0_complete_QC_review_pending')
                self.assertTrue(state['pilot_passed']);combined_mock.assert_called_once()
                self.assertEqual(publishing.call_count,4)
                self.assertEqual(len(json.loads((record/'regenerated-ugr-bold-overrides.json').read_text())),249)
        self.assertEqual(sentinel.read_bytes(),b'previously completed Trust output')

    def test_pilot_pass_expands_automatically_and_publishes_complete_combined_registry(self):
        with tempfile.TemporaryDirectory() as d:self.workflow(d)

    def test_failed_pilot_prevents_production_expansion(self):
        with tempfile.TemporaryDirectory() as d:self.workflow(d,fail_pilot=True)


if __name__=='__main__':unittest.main()
