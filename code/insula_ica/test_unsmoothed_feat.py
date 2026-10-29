from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

import run_unsmoothed_feat as s


class UnsmoothedTests(unittest.TestCase):
    def fixture(self,root,present=False):
        feat=root/'old.feat';feat.mkdir();(root/'fsf').mkdir();record=root/'record';(record/'qc').mkdir(parents=True)
        conf=root/'confounds.txt';np.savetxt(conf,np.zeros((6,1)));(feat/'confoundevs.txt').write_text(conf.read_text())
        for name,matrix in [('design.mat',np.zeros((6,2))),('design.con',np.ones((1,2)))]:
            (feat/name).write_text('/Matrix\n'+'\n'.join(' '.join(str(x) for x in row) for row in matrix)+'\n')
        bold=root/'preproc_bold.nii.gz';event=root/'task_ev.txt';event.write_text('0 1 1\n')
        nib.save(nib.Nifti1Image(np.ones((5,5,5),dtype='float32'),np.eye(4)),feat/'mean_func.nii.gz')
        if present:nib.save(nib.Nifti1Image(np.ones((5,5,5,6),dtype='float32'),np.eye(4)),bold)
        text=f'set fmri(smooth) 5\nset fmri(outputdir) "{feat}"\nset feat_files(1) "{bold}"\nset confoundev_files(1) "{conf}"\nset fmri(evs_orig) 2\nset fmri(shape1) 3\nset fmri(custom1) "{event}"\nset fmri(shape2) 10\nset fmri(temphp_yn) 0\nset fmri(prewhiten_yn) 1\n'
        fsf=feat/'design.fsf';fsf.write_text(text)
        row=dict(key='test',task='trust',active_feat=str(feat),source_fsf=str(fsf),confounds=str(conf),confounds_sha256=s.g.core.sha(conf),nvol=6)
        return row,text,record

    def test_smoothing_change_preserves_all_other_model_settings(self):
        text='set fmri(smooth) 5\nset fmri(outputdir) "/old"\nset feat_files(1) "/bold"\nset confoundev_files(1) "/conf"\nset fmri(temphp_yn) 0\nset fmri(prewhiten_yn) 1\n'
        new=s.unsmoothed_fsf(text,'/bold','/conf','/new')
        self.assertEqual(new.replace('fmri(smooth) 0','fmri(smooth) 5').replace('/new','/old'),text)

    def test_missing_source_bold_does_not_fall_back_to_smoothed_filtered_data(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row,text,_=self.fixture(root)
            (Path(row['active_feat'])/'filtered_func_data.nii.gz').touch()
            with patch.object(s.g.residuals,'validate_model'),patch.object(s.subprocess,'run') as subprocess:
                plan=s.plan(row,root)
                subprocess.assert_not_called()
            self.assertEqual(plan['status'],'missing_source_inputs')
            self.assertEqual(plan['missing_inputs'],[str(root/'preproc_bold.nii.gz')])
            self.assertEqual(Path(row['source_fsf']).read_text(),text)

    def test_ready_plan_uses_exact_saved_bold_and_nuisance_paths(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row,text,_=self.fixture(root,True)
            with patch.object(s.g.residuals,'validate_model'):plan=s.plan(row,root)
            self.assertEqual(plan['status'],'ready');self.assertEqual(plan['source_bold'],str(root/'preproc_bold.nii.gz'))
            settings=s.g.core.fsl_settings(plan['unsmoothed_fsf'])
            self.assertEqual(settings['fmri(smooth)'],'0');self.assertEqual(settings['fmri(temphp_yn)'],'0')
            self.assertEqual(settings['fmri(prewhiten_yn)'],'1');self.assertEqual(settings['confoundev_files(1)'],row['confounds'])
            self.assertEqual(Path(row['source_fsf']).read_text(),text)

    def test_partial_feat_is_preserved_without_executing_the_wrapper(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row,_,record=self.fixture(root,True)
            with patch.object(s.g.residuals,'validate_model'):plan=s.plan(row,root)
            Path(plan['target']).mkdir(parents=True)
            with patch.object(s.subprocess,'run') as execute:
                with self.assertRaisesRegex(ValueError,'Partial FEAT exists'):s.fit(plan,root,record)
                execute.assert_not_called()

    def test_completed_sm0_fit_is_reused_with_unchanged_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row,_,record=self.fixture(root,True)
            with patch.object(s.g.residuals,'validate_model'):plan=s.plan(row,root)
            target=Path(plan['target']);target.mkdir(parents=True);output=target/'residual';output.write_bytes(b'frozen output')
            identity=s.g.core.identity(output,False)
            signature=dict(bold=plan['source_bold_identity'],events=plan['event_identities'],original_fsf_sha256=plan['source_fsf_sha256'],
                confounds_sha256=plan['confounds_sha256'],unsmoothed_fsf_sha256=s.g.core.sha(plan['unsmoothed_fsf']))
            (target/'sm0_complete.json').write_text(json.dumps(dict(signature=signature,outputs=[identity],qc=dict(key='test'))))
            with patch.object(s.subprocess,'run') as execute:
                result=s.fit(plan,root,record);execute.assert_not_called()
            self.assertTrue(result['reused_completed_feat']);self.assertEqual(s.g.core.identity(output,False),identity)


if __name__=='__main__':unittest.main()
