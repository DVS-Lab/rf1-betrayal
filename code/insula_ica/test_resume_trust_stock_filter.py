from contextlib import ExitStack
import copy
import json
from pathlib import Path
import signal
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

import resume_trust_stock_filter as r
from diagnose_trust_registration import sbref_selection


class ResumeTests(unittest.TestCase):
    def test_handoff_signals_only_scheduler_and_retains_running_feat(self):
        with tempfile.TemporaryDirectory() as directory,ExitStack() as p:
            root=Path(directory);record=root/'record';record.mkdir();output=root/'output'
            parent=dict(pid=100,ppid=1,start='10',state='S',cwd=str(root),args=['/python','rerun_trust_references.py'])
            child=dict(pid=101,ppid=100,start='11',state='R',cwd=str(root),
                       args=['bash','/fsl/feat',str(output/'sub-1/run.fsf')])
            table={100:parent,101:child};signals=[];pending=False
            def send(pid,sig):
                nonlocal pending
                signals.append((pid,sig));self.assertEqual(pid,100)
                if sig==signal.SIGSTOP:parent['state']='T'
                if sig==signal.SIGTERM:pending=True
                if sig==signal.SIGCONT and pending:table.pop(100)
            p.enter_context(patch.object(r,'ROOT',root))
            p.enter_context(patch.object(r,'processes',side_effect=lambda:copy.deepcopy(table)))
            p.enter_context(patch.object(r,'stopped_threads',return_value=True))
            p.enter_context(patch.object(r.os,'kill',side_effect=send))
            inherited=r.handoff(output,record)
            self.assertEqual([p['pid'] for p in inherited],[101])
            self.assertIn(101,table)
            self.assertEqual(signals,[(100,signal.SIGSTOP),(100,signal.SIGTERM),(100,signal.SIGCONT)])

    def test_active_fmriprep_is_allowed_to_finish_before_handoff(self):
        with tempfile.TemporaryDirectory() as directory,ExitStack() as p:
            root=Path(directory);record=root/'record';record.mkdir();output=root/'output'
            parent=dict(pid=100,ppid=1,start='10',state='S',cwd=str(root),args=['/python','rerun_trust_references.py'])
            child=dict(pid=101,ppid=100,start='11',state='R',cwd=str(root),args=['apptainer','run','image.sif'])
            table={100:parent,101:child};signals=[];pending=False
            def send(pid,sig):
                nonlocal pending
                self.assertEqual(pid,100);signals.append(sig)
                if sig==signal.SIGSTOP:parent['state']='T'
                if sig==signal.SIGTERM:pending=True
                if sig==signal.SIGCONT:
                    if pending:table.pop(100)
                    else:
                        parent['state']='S'
                        child['args']=['bash','/fsl/feat',str(output/'sub-1/run.fsf')]
            p.enter_context(patch.object(r,'ROOT',root))
            p.enter_context(patch.object(r,'processes',side_effect=lambda:copy.deepcopy(table)))
            p.enter_context(patch.object(r,'stopped_threads',return_value=True))
            p.enter_context(patch.object(r.os,'kill',side_effect=send))
            p.enter_context(patch.object(r.time,'sleep'))
            r.handoff(output,record)
            self.assertEqual(signals[:2],[signal.SIGSTOP,signal.SIGCONT])
            self.assertEqual(signals.count(signal.SIGTERM),1)

    def test_pid_reuse_and_unrelated_children_are_rejected(self):
        identity=dict(pid=1,start='old')
        self.assertFalse(r.alive(identity,{1:dict(pid=1,start='new',state='R')}))
        table={2:dict(pid=2,ppid=1,state='R',args=['feat','/other/run.fsf'])}
        self.assertIsNone(r.feat_children(table,identity,Path('/expected')))

    def test_two_preprocessing_slots_reserve_ram_and_release_it(self):
        with patch.object(r.MemoryBudget,'available',return_value=100):
            budget=r.MemoryBudget(2);budget.acquire('fmriprep');budget.acquire('fmriprep')
            self.assertEqual((budget.reserved,budget.fmriprep),(56,2))
            budget.release('fmriprep');budget.release('fmriprep');budget.acquire('feat_qc')
            self.assertEqual((budget.reserved,budget.fmriprep),(4,0))
            budget.release('feat_qc');self.assertEqual(budget.reserved,0)

    def test_actual_archived_stock_logs_verify_with_scoped_input_view(self):
        archive=Path(__file__).parent/'records/linux2-trust-progress-20261028T195956Z-3cc48c'
        for subject in ('sub-10402','sub-10478','sub-10596'):
            with self.subTest(subject=subject),tempfile.TemporaryDirectory() as directory:
                case=Path(directory);plan=json.loads((archive/subject/'case.json').read_text())
                row=dict(plan['original_row']);key=row['key']
                row['new_bold']=str(case/'fmriprep-24'/subject/'ses-01/func'/(key+'_bold.nii.gz'))
                view=case/'bids-view'/subject/'ses-01/func';view.mkdir(parents=True)
                for name in row['raw_images']:(view/Path(name).name).touch()
                refs=json.loads((archive/subject/'reference-files.json').read_text())
                for ref in refs:(view/Path(ref['view']).name).touch()
                (case/(subject+'.fmriprep.log')).write_text((archive/subject/'reference-log-excerpt.txt').read_text())
                result=sbref_selection([row])[0]
                self.assertEqual(result['status'],'same_run')
                self.assertEqual(result['evidence_type'],'explicit_selection_log_with_verified_single_run_input_view')
                # A view containing another BOLD run destroys unambiguous attribution.
                (view/(key.replace('_run-2','_run-1')+'_bold.nii.gz')).touch()
                self.assertEqual(sbref_selection([row])[0]['status'],'unknown_no_retained_log_evidence')

    def test_wrong_run_selection_is_never_accepted_by_new_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory);subject='sub-10001';key=subject+'_ses-01_task-trust_run-2'
            view=case/'bids-view'/subject/'ses-01/func';view.mkdir(parents=True)
            name=key+'_echo-1_part-mag_bold.nii.gz';(view/name).touch()
            other=key.replace('_run-2','_run-1')+'_echo-1_part-mag_sbref.nii.gz';(view/other).touch()
            (case/(subject+'.fmriprep.log')).write_text('Using single-band reference file(s) '+other+'.\n')
            row=dict(subject=subject,key=key,raw_images={name:{}},
                new_bold=str(case/'fmriprep-24'/subject/'ses-01/func'/'bold.nii.gz'))
            self.assertEqual(sbref_selection([row])[0]['status'],'unknown_no_retained_log_evidence')

    def test_completed_qc_is_reused_without_fit_or_residual_array_read(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source';record=root/'record';case=root/'case'
            for folder in (source/'qc',record/'qc',case/'feat_reconstruction/run.feat/stats',root/'historical'):
                folder.mkdir(parents=True)
            target=case/'feat_reconstruction/run.feat';old=root/'historical'
            bold=case/'bold.nii.gz';image=nib.Nifti1Image(np.zeros((2,2,2,3),np.float32),np.eye(4))
            nib.save(image,bold);nib.save(image,target/'stats/res4d.nii.gz')
            nib.save(nib.Nifti1Image(np.ones((2,2,2),np.float32),np.eye(4)),target/'mask.nii.gz')
            conf=root/'preserved.txt';conf.write_text('1\n2\n3\n')
            for folder in (old,target):
                (folder/'confoundevs.txt').write_bytes(conf.read_bytes())
                for name in ('design.mat','design.con'):
                    (folder/name).write_text('/NumWaves 1\n/NumPoints 1\n/Matrix\n1\n')
            (target/'historical_comparison.json').write_text(json.dumps(dict(ncontrasts=1,
                missing_original=[],missing_reconstructed=[],comparisons=[dict(grid_matches=True)]*2)))
            qc=dict(finite_nonconstant_residuals=True,residual=str(target/'stats/res4d.nii.gz'),review_flags=[])
            (source/'qc/run.json').write_text(json.dumps(qc));(source/'qc/run.png').write_bytes(b'previous montage')
            row=dict(key='run',new_bold=str(bold),original_feat=str(old),confounds=str(conf),confounds_sha256=r.original.core.sha(conf))
            with patch.object(r.original.rebuild,'verify_outputs') as verify:
                self.assertEqual(r.reuse_qc(row,case,source,record),qc)
                verify.assert_called_once_with([row])
            self.assertEqual((record/'qc/run.json').read_bytes(),(source/'qc/run.json').read_bytes())
            # Touching a fitted input invalidates the old QC cache.
            (target/'confoundevs.txt').write_text('1\n2\n3\n')
            self.assertIsNone(r.reuse_qc(row,case,source,record))

    def test_resume_pipeline_skips_completed_fmriprep_and_fit(self):
        with tempfile.TemporaryDirectory() as directory,ExitStack() as p:
            root=Path(directory);source=root/'source';source.mkdir();record=root/'new-record';record.mkdir()
            subject='sub-10001';key=subject+'_ses-01_task-trust_run-2'
            output=root/'derivatives/insula_ica/trust-stock-filter-fixture';case=output/subject;case.mkdir(parents=True)
            registry=root/'code/insula_ica/records'/r.original.SOURCE;registry.mkdir(parents=True)
            conf=root/'confounds.txt';conf.write_text('1\n2\n')
            row=dict(subject=subject,key=key,run=2,confounds=str(conf),confounds_sha256=r.original.core.sha(conf),
                     new_bold='/previous/fmriprep-24/sub-10001/ses-01/func/example_bold.nii.gz')
            (registry/'container.json').write_text(json.dumps(dict(sha256='same-image')))
            (registry/'input-registry.json').write_text(json.dumps([row]))
            (source/'state.json').write_text(json.dumps(dict(output_root=str(output))))
            (source/'reference-selection-audit.json').write_text('[]')
            (source/'original-fmriprep.sh').write_text('# unchanged original script\n')
            stock=root/'stock.json';stock.write_text('{}')
            (source/'fmriprep_config.json').write_text('{}');(case/'fmriprep_config.json').write_text('{}')
            (source/'cases').mkdir()
            (source/'cases'/(subject+'.json')).write_text(json.dumps(dict(argv=['previous-command'],anatomical_inputs=[],state='preflight')))
            (case/(subject+'.fmriprep.log')).write_text('ANAT Stage 4: Found pre-computed registrations\nfMRIPrep finished successfully!\n')
            (case/'anatomy-preflight.json').write_text(json.dumps(dict(anatomical_registration_nodes_absent=True)))
            repair=root/'code/insula_ica/records/linux2-registration-fix-sub-10834-20261028T151154Z-acdf91/review-decision.json'
            repair.parent.mkdir(parents=True);repair.write_text(json.dumps(dict(selected_inputs=[])))
            p.enter_context(patch.object(r,'ROOT',root));p.enter_context(patch.object(r.original,'STOCK',stock))
            p.enter_context(patch.object(r.original.rebuild,'sha',return_value='same-image'))
            p.enter_context(patch.object(r.shutil,'which',side_effect=lambda name:'/bin/'+name))
            p.enter_context(patch.dict(r.os.environ,FSLDIR='/fixture'))
            p.enter_context(patch.object(r.original,'choose_runs',return_value=[row]))
            p.enter_context(patch.object(r,'handoff',return_value=[]));p.enter_context(patch.object(r,'case_processes',return_value=[]))
            p.enter_context(patch.object(r,'sbref_selection',return_value=[dict(status='same_run')]))
            p.enter_context(patch.object(r,'reuse_qc',return_value=dict(review_flags=[])))
            execute=p.enter_context(patch.object(r.original,'execute'))
            fit=p.enter_context(patch.object(r.original.rebuild,'fit_run'))
            qc=p.enter_context(patch.object(r,'qc_run'))
            state=dict(subjects={});r.run(SimpleNamespace(record=source,jobs=16,fmriprep_jobs=2),record,state)
            self.assertEqual(state['state'],'functional_reruns_complete_qc_pending')
            self.assertEqual(state['completed'],1)
            self.assertTrue(state['subjects'][subject]['fmriprep_reused'])
            self.assertTrue(state['subjects'][subject]['numerical_qc_reused'])
            execute.assert_not_called();fit.assert_not_called();qc.assert_not_called()


if __name__=='__main__':
    unittest.main()
