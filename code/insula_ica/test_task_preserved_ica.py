import unittest
from contextlib import ExitStack,nullcontext
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch
import nibabel as nib
import numpy as np
import task_preserved_ica as s


class TaskPreservedTests(unittest.TestCase):
    def test_nuisance_only_preserves_conditional_task_signal(self):
        rng=np.random.default_rng(7);q,_=np.linalg.qr(rng.normal(size=(100,8)))
        q-=q.mean(axis=0);q,_=np.linalg.qr(q)
        nuisance=q[:,:3];task=q[:,3:5]+nuisance[:,:2]*.7
        design=np.column_stack([task,nuisance]);betas=rng.normal(size=(4,5))
        noise=rng.normal(size=(4,3))@q[:,5:8].T
        y=betas@design.T+noise+10000
        results=s.ols_pair(y,design,2)
        np.testing.assert_allclose(results['full_model_ols'],noise,atol=1e-10)
        np.testing.assert_allclose(results['confounds_only_ols'],noise+betas[:,:2]@q[:,3:5].T,atol=1e-10)
        self.assertGreater(np.linalg.norm(results['confounds_only_ols']-results['full_model_ols']),1)
        np.testing.assert_allclose(results['confounds_only_ols']@nuisance,0,atol=1e-10)

    def test_mean_offset_does_not_change_the_matched_regressions(self):
        rng=np.random.default_rng(19);x=rng.normal(size=(100,9));x-=x.mean(axis=0)
        y=rng.normal(size=(20,100))
        a,b=s.ols_pair(y,x,4),s.ols_pair(y+500,x,4)
        for mode in s.MODES:np.testing.assert_allclose(a[mode],b[mode],atol=1e-10)

    def test_no_task_columns_or_missing_nuisance_block_is_rejected(self):
        with self.assertRaises(ValueError):s.ols_pair(np.ones((3,100)),np.ones((100,5)),5)
        with self.assertRaises(ValueError):s.ols_pair(np.ones((3,99)),np.ones((100,5)),2)

    def test_full_cohort_produces_six_matched_fits_without_upstream_jobs(self):
        with tempfile.TemporaryDirectory() as d,ExitStack() as p:
            base=Path(d).resolve();root=base/'analysis';root.mkdir();(root/'run_arrays').mkdir();record=base/'record';record.mkdir()
            for name in ('runs','errors','attempts','figures'):(record/name).mkdir()
            source=base/'source';source.mkdir();parent=base/'existing';parent.mkdir()
            (source/'state.json').write_text(json.dumps(dict(output_root=str(parent))))
            original=Path(__file__).parent/'records/linux2-paired-group-ica-20261028T215022Z-679a3e/accepted-inputs.json'
            rows=json.loads(original.read_text())
            for row in rows:row['nvol']=20
            (source/'accepted-inputs.json').write_text(json.dumps(rows))
            reference=nib.Nifti1Image(np.ones((5,5,5),dtype='float32'),np.eye(4));masks={roi:np.ones(reference.shape,bool) for roi in s.g.core.ROIS}
            p.enter_context(patch.object(s.resume,'verify_prepared',return_value=({},source,{},masks,None,reference)))
            p.enter_context(patch.object(s.g.core,'software',return_value=dict(fsl='fixture')))
            p.enter_context(patch.object(s.g.inputs.MemoryGate,'reserve',side_effect=lambda amount:nullcontext()))
            p.enter_context(patch.object(s.g.inputs,'ensure_native',side_effect=AssertionError('No FILM')))
            def prepare(row,root,record,union,ref):
                rng=np.random.default_rng(row['run'])
                for mode in s.MODES:np.save(root/'run_arrays'/(row['key']+'_'+mode+'.npy'),rng.normal(size=(125,20)).astype('float32'))
                (root/'run_arrays'/(row['key']+'_complete.json')).write_text('{}')
            p.enter_context(patch.object(s,'prepare',side_effect=prepare))
            def fit(plan,root,masks,slices,ref,record,**kwargs):
                mode,roi,dim=plan;folder=root/'fits'/mode/roi;folder.mkdir(parents=True)
                path=folder/'maps.nii.gz';s.g.core.save_image(path,np.random.default_rng(5).normal(size=(5,5,5,5)).astype('float32'),ref)
                self.assertTrue(kwargs['variance_normalization'])
                return dict(state='complete',oriented_maps=str(path),mode=mode,roi=roi,dimensionality=5)
            fits=p.enter_context(patch.object(s.resume,'retry_fit',side_effect=fit))
            p.enter_context(patch.object(s.g,'component_figure',side_effect=lambda maps,mask,bg,aff,path,label:path.write_bytes(b'fixture')))
            state=dict(runs={},fits={});s.run(SimpleNamespace(source_record=source,jobs=64),root,record,state)
            self.assertEqual(state['state'],'task_preserved_sensitivity_complete_review_pending')
            self.assertEqual(len(state['runs']),490);self.assertEqual(fits.call_count,6)
            self.assertEqual(len(list((record/'figures').glob('*.png'))),6)
            self.assertEqual(len(json.loads((record/'task_removal_spatial_comparisons.json').read_text())),3)


if __name__=='__main__':unittest.main()
