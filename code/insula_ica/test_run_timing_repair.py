from contextlib import ExitStack
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import run_timing as timing
import run_unsmoothed_feat as sm0
import test_unsmoothed_feat as fixture
from test_insula_ica import vest
import repair_ugr_timing as repair


class TimingTests(unittest.TestCase):
    def test_raw_consensus_checks_all_echoes_and_units_without_altering_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);key='sub-1_ses-01_task-ugr_run-1'
            for echo in range(1,5):
                im=nib.Nifti1Image(np.zeros((2,2,2,6),dtype='float32'),np.eye(4))
                im.header.set_xyzt_units('mm','msec');im.header.set_zooms((1,1,1,1710))
                stem=key+f'_echo-{echo}_part-mag_bold'
                nib.save(im,folder/(stem+'.nii.gz'))
                (folder/(stem+'.json')).write_text(json.dumps(dict(RepetitionTime=1.71,SliceTiming=[0,.5,1.])))
            before={str(p):p.read_bytes() for p in folder.iterdir()}
            self.assertEqual(timing.acquisition_timing(folder,key)['tr_seconds'],1.71)
            self.assertEqual(before,{str(p):p.read_bytes() for p in folder.iterdir()})
            side=folder/(key+'_echo-4_part-mag_bold.json');side.write_text(json.dumps(dict(RepetitionTime=1.615)))
            with self.assertRaisesRegex(ValueError,'JSON TR'):timing.acquisition_timing(folder,key)

    def test_FSF_stamping_changes_only_the_run_TR_and_refuses_ambiguity(self):
        text='set fmri(tr) 1.615000\nset fmri(smooth) 5\nset fmri(custom1) "/events_seconds.txt"\n'
        self.assertEqual(timing.stamp(text,1.71).replace('1.71','1.615000'),text)
        with self.assertRaisesRegex(ValueError,'exactly one'):timing.stamp(text+text,1.71)

    def test_both_project_launchers_validate_before_skipping_and_stamp_both_models(self):
        for task in ('ugr','trust'):
            text=(Path(__file__).resolve().parents[1]/f'L1stats-{task}.sh').read_text()
            self.assertLess(text.index('run_timing.py" resolve'),text.index('cluster_mask_zstat1'))
            self.assertIn('check-existing --fsf',text)
            self.assertEqual(text.count('stamp --fsf'),2)

    def simulate_corrected_fit(self,root,header_only=False):
        row,text,record=fixture.UnsmoothedTests().fixture(root)
        (root/'logs').mkdir()
        fsf=Path(row['source_fsf']);fsf.write_text(text+'set fmri(tr) 1.615\nset fmri(evs_real) 1\n')
        row.update(task='ugr',tr=1.615)
        correction=dict(repair.CORRECTION,key=row['key'])
        bold=root/'new_1p71.nii.gz';image=nib.Nifti1Image(np.ones((5,5,5,6),dtype='float32'),np.eye(4));image.header.set_zooms((1,1,1,1.71));nib.save(image,bold)
        override=dict(path=str(bold),identity=sm0.g.core.identity(bold,False),timing_correction=correction)
        with patch.object(sm0.g.residuals,'validate_model'):item=sm0.plan(row,root,override)
        self.assertEqual(float(sm0.g.core.fsl_settings(item['unsmoothed_fsf'])['fmri(tr)']),1.71)
        original_bytes=fsf.read_bytes();design=np.zeros((6,2))
        if not header_only:design[:,0]=np.linspace(-1,1,6)
        def execute(argv,**kwargs):
            if argv[0]=='feat_model':
                base=Path(argv[1]);vest(base.with_suffix('.mat'),design);vest(base.with_suffix('.con'),np.ones((1,2)))
            else:
                out=Path(item['target']);(out/'stats').mkdir(parents=True)
                (out/'design.fsf').write_text(Path(item['unsmoothed_fsf']).read_text())
                vest(out/'design.mat',design);vest(out/'design.con',np.ones((1,2)))
                (out/'confoundevs.txt').write_text(Path(row['confounds']).read_text());(out/'stats/dof').write_text('50')
                values=np.broadcast_to(np.arange(6,dtype='float32'),(5,5,5,6)).copy()
                im=nib.Nifti1Image(values,np.eye(4));im.header.set_zooms((1,1,1,1.71))
                for path in ('stats/res4d.nii.gz','filtered_func_data.nii.gz'):nib.save(im,out/path)
                for path in ('mask.nii.gz','mean_func.nii.gz','stats/threshac1.nii.gz','stats/sigmasquareds.nii.gz','stats/pe1.nii.gz','stats/pe2.nii.gz'):
                    nib.save(nib.Nifti1Image(np.ones((5,5,5),dtype='float32'),np.eye(4)),out/path)
        masks={roi:np.ones((5,5,5),bool) for roi in ('left','right','bilateral')}
        with patch.object(sm0.subprocess,'run',side_effect=execute),patch.object(sm0.g.core,'anatomical_masks',return_value=(masks,{},'insula')),patch.dict(os.environ,FSLDIR='/fixture'):
            if header_only:
                with self.assertRaisesRegex(ValueError,'header-only'):sm0.fit(item,root,record)
                self.assertFalse((Path(item['target'])/'sm0_complete.json').exists())
            else:
                result=sm0.fit(item,root,record);self.assertEqual(result['state'],'complete')
                qc=json.loads((record/'qc'/(row['key']+'.json')).read_text())
                self.assertFalse(qc['complete_design_preserved']);self.assertTrue(qc['task_design_regenerated_for_TR'])
                self.assertTrue(qc['exact_original_nuisance_matrix']);self.assertTrue(qc['contrast_definitions_preserved'])
        self.assertEqual(fsf.read_bytes(),original_bytes)

    def test_corrected_design_is_validated_against_regenerated_model_not_historical_model(self):
        with tempfile.TemporaryDirectory() as d:self.simulate_corrected_fit(Path(d))

    def test_header_only_TR_repair_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:self.simulate_corrected_fit(Path(d),True)


if __name__=='__main__':unittest.main()
