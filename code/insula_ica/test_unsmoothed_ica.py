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

import run_unsmoothed_ica as u
import run_group_ica as g


class UnsmoothedICATests(unittest.TestCase):
    def block_rows(self):
        return [dict(subject=str(10000+j),task=task,run=1,nvol=12,
                     tr=1.710 if j==0 and task=='ugr' else 1.615,
                     key=f'sub-{10000+j}_ses-01_task-{task}_run-1')
                for j in range(132) for task in ('trust','ugr')]

    def masks(self):
        left=np.zeros((3,3,3),bool);left[0]=True;right=np.zeros_like(left);right[2]=True
        return dict(left=left,right=right,bilateral=left|right)

    def test_physical_4D_inputs_preserve_both_TRs_and_specific_ROI_masks(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);record=root/'record';record.mkdir();(root/'masks').mkdir();(root/'run_arrays').mkdir()
            masks=self.masks();union=masks['bilateral'];rows=self.block_rows();rng=np.random.default_rng(7)
            reference=nib.Nifti1Image(np.ones((3,3,3),dtype='float32'),np.eye(4))
            for row in rows:
                for mode in g.inputs.MODES:np.save(root/'run_arrays'/(row['key']+'_'+mode+'.npy'),rng.normal(size=(int(union.sum()),12)).astype('float32'))
                (root/'run_arrays'/(row['key']+'_complete.json')).write_text('{}')
            g.sessions(rows,root,record,masks,union,reference,amplitude=1000,preserve_timing=True)
            for mode in g.inputs.MODES:
                for roi in g.core.ROIS:
                    folder=root/'sessions'/mode/roi
                    self.assertEqual(len((folder/'inputs.txt').read_text().splitlines()),264)
                    im=nib.load(folder/'sub-10000_task-ugr.nii.gz');other=nib.load(folder/'sub-10000_task-trust.nii.gz')
                    self.assertAlmostEqual(float(im.header.get_zooms()[3]),1.710,places=5)
                    self.assertAlmostEqual(float(other.header.get_zooms()[3]),1.615,places=5)
                    self.assertEqual(np.count_nonzero(np.asarray(im.dataobj)[~masks[roi]]),0)
                    self.assertTrue(np.isfinite(np.asarray(im.dataobj)).all())

    def test_mixed_TRs_within_a_block_are_not_silently_relabelled(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'masks').mkdir();record=root/'record';record.mkdir()
            rows=self.block_rows();rows.append(dict(rows[0],run=2,tr=1.710,key='mixed_run_2'))
            masks=self.masks();reference=nib.Nifti1Image(np.ones((3,3,3),dtype='float32'),np.eye(4))
            with self.assertRaisesRegex(ValueError,'Different TRs within'):
                g.sessions(rows,root,record,masks,masks['bilateral'],reference,preserve_timing=True)

    def test_command_keeps_concat_and_VN_without_common_TR_or_frequency_report(self):
        argv=g.melodic_command('joint_bold','left',5,Path('/new'),output_folder=Path('/new/retry'),
            maxit=10000,variance_normalization=True,preserve_timing=True)
        self.assertEqual(argv[argv.index('-a')+1],'concat');self.assertNotIn('--vn',argv)
        self.assertNotIn('--tr=1.615',argv);self.assertNotIn('--report',argv)
        self.assertIn('--disableMigp',argv);self.assertIn('--update_mask',argv);self.assertIn('--nl=pow3',argv)

    def test_source_registry_uses_all_490_completed_fits_and_requires_corrected_design_QC(self):
        original=json.loads((g.core.ROOT/'code/insula_ica/records/linux2-paired-group-ica-20261028T215022Z-679a3e/accepted-inputs.json').read_text())
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'source';source.mkdir();(source/'qc').mkdir();base=root/'base';base.mkdir()
            (base/'accepted-inputs.json').write_text(json.dumps(original));featroot=root/'sm0';featroot.mkdir()
            (root/'code').mkdir();(root/'code/sublist_n132.txt').write_text('\n'.join(sorted({r['subject'] for r in original})))
            states={};corrected_qc=None;corrected_marker=None
            for row in original:
                key=row['key'];fixed=key=='sub-10817_ses-01_task-ugr_run-1'
                feat=featroot/'feat'/(key+'_sm-0.feat');feat.mkdir(parents=True)
                fsf=feat/'design.fsf';fsf.write_text(f'set fmri(smooth) 0\nset fmri(tr) {1.710 if fixed else 1.615}\n')
                qc=dict(finite_residuals=True,exact_original_nuisance_matrix=True,complete_design_preserved=not fixed)
                if fixed:qc['task_design_regenerated_for_TR']=True
                path=source/'qc'/(key+'.json');path.write_text(json.dumps(qc))
                marker=feat/'sm0_complete.json';marker.write_text(json.dumps(dict(qc=qc,outputs=[g.core.identity(fsf)])))
                if fixed:corrected_qc,corrected_marker=path,marker
                states[key]=dict(state='complete')
            (source/'state.json').write_text(json.dumps(dict(state='sm0_FEAT_complete_QC_pending',runs=states,source_record=str(base),output_root=str(featroot))))
            with patch.object(u,'ROOT',root):
                rows=u.select_inputs(source)
                self.assertEqual(len(rows),490);fixed=next(r for r in rows if r['key']=='sub-10817_ses-01_task-ugr_run-1')
                self.assertEqual(fixed['tr'],1.710);self.assertEqual(fixed['original_feat'],fixed['active_feat'])
                self.assertNotEqual(fixed['manuscript_original_feat'],fixed['active_feat'])
                qc=json.loads(corrected_qc.read_text());qc['task_design_regenerated_for_TR']=False;corrected_qc.write_text(json.dumps(qc))
                marker=json.loads(corrected_marker.read_text());marker['qc']=qc;corrected_marker.write_text(json.dumps(marker))
                with self.assertRaisesRegex(ValueError,'task design not verified'):u.select_inputs(source)

    def orchestration(self,root,fail_fixed):
        record=root/'record';record.mkdir();analysis=root/'analysis';analysis.mkdir()
        for name in ('runs','figures','attempts','errors'):(record/name).mkdir()
        old=root/'old';old.mkdir();native=old/'res4d.nii.gz';native.touch()
        reference=nib.Nifti1Image(np.ones((3,3,3),dtype='float32'),np.eye(4));nib.save(reference,old/'mean_func.nii.gz')
        rows=[dict(self.block_rows()[j%264],key=f'run-{j}',active_feat=str(old),native_residual=str(native)) for j in range(490)]
        atlas=root/'fsl/data/atlases/HarvardOxford';atlas.mkdir(parents=True)
        (atlas/'HarvardOxford-cort-maxprob-thr25-2mm.nii.gz').touch();(atlas.parent/'HarvardOxford-Cortical.xml').touch()
        state=dict(runs={},fits={});calls=[]
        def fit(plan,*a,**kw):
            calls.append((plan,kw));self.assertTrue(kw['variance_normalization']);self.assertTrue(kw['preserve_timing'])
            mode,roi,dim=plan
            if fail_fixed and plan==('joint_bold','left',7):raise ValueError('No convergence after five deterministic starts')
            self.assertEqual(kw['threads'],8 if dim==0 else 4)
            return dict(mode=mode,roi=roi,dimensionality=dim,state='complete')
        with ExitStack() as p,redirect_stdout(io.StringIO()):
            p.enter_context(patch.dict(os.environ,FSLDIR=str(root/'fsl')))
            p.enter_context(patch.object(u,'select_inputs',return_value=rows))
            p.enter_context(patch.object(g.core,'software',return_value={}))
            p.enter_context(patch.object(g.core,'anatomical_masks',return_value=(self.masks(),{},'insula')))
            p.enter_context(patch.object(g.inputs,'inspect_source'))
            p.enter_context(patch.object(g.shutil,'disk_usage',return_value=SimpleNamespace(free=2**40)))
            p.enter_context(patch.object(g.residuals,'residual_qc',return_value=dict(scale='voxelwise_prewhitened',max_whitened_decomposition_relative_error=0)))
            p.enter_context(patch.object(g.inputs.MemoryGate,'reserve',side_effect=lambda amount:nullcontext()))
            p.enter_context(patch.object(g.inputs,'prepare_run',return_value=dict(native_recovered_by_film=False)))
            p.enter_context(patch.object(g,'fixed_masks',return_value=self.masks()))
            sessions=p.enter_context(patch.object(g,'sessions',return_value=((slice(0,3),)*3,np.eye(4))))
            p.enter_context(patch.object(u.retry,'retry_fit',side_effect=fit))
            review=p.enter_context(patch.object(g,'review_outputs'))
            partial=p.enter_context(patch.object(u,'render_partial'))
            u.run(SimpleNamespace(source_record=Path('/validated'),jobs=64),analysis,record,state)
            self.assertEqual(len(calls),24);self.assertTrue(all(dim!=0 for (m,r,dim),kw in calls[:18]))
            self.assertTrue(all(dim==0 for (m,r,dim),kw in calls[18:]))
            self.assertEqual(len(state['runs']),490)
            self.assertTrue(sessions.call_args.kwargs['preserve_timing']);self.assertEqual(sessions.call_args.kwargs['amplitude'],1000)
            if fail_fixed:
                self.assertEqual(state['state'],'unsmoothed_ICA_incomplete_review_errors');partial.assert_called_once();review.assert_not_called()
            else:self.assertEqual(state['state'],'unsmoothed_ICA_complete_review_pending');review.assert_called_once()

    def test_orchestration_runs_all_24_and_uses_prepared_native_sources(self):
        with tempfile.TemporaryDirectory() as d:self.orchestration(Path(d),False)

    def test_automatic_fits_still_run_after_fixed_nonconvergence(self):
        with tempfile.TemporaryDirectory() as d:self.orchestration(Path(d),True)


if __name__=='__main__':unittest.main()
