import json
from pathlib import Path
import tempfile
import unittest

import nibabel as nib
import numpy as np
from scipy import stats

import explore_completed_ica as e
import run_group_ica as g
from test_insula_ica import vest


class ExplorationTests(unittest.TestCase):
    def test_joint_wholebrain_projection_is_exact_existing_GLS_residual_not_new_OLS(self):
        rng=np.random.default_rng(7);t=80
        data=rng.normal(size=(27,t))+1000;design=rng.normal(size=(t,8));beta=rng.normal(size=(27,8))
        operator=rng.normal(size=(5,t));residual=(data-data.mean(axis=1,keepdims=True))-beta@design.T
        np.testing.assert_allclose(e.joint_coefficients(operator,data,design,beta),operator@residual.T,atol=1e-10)
        different=(data-data.mean(axis=1,keepdims=True))@np.linalg.pinv(design).T
        self.assertGreater(np.linalg.norm(beta-different),1)

    def test_spatial_projection_fits_all_maps_and_removes_only_spatial_intercept(self):
        rng=np.random.default_rng(8);maps=rng.normal(size=(90,5));maps[:,1]+=.7*maps[:,0]
        tc=rng.normal(size=(5,70));offset=rng.normal(size=70)+1000
        weights=e.weights(maps)
        np.testing.assert_allclose(weights@(maps@tc+offset),tc,atol=1e-11)
        np.testing.assert_allclose(weights.sum(axis=1),0,atol=1e-12)

    def test_temporal_projection_is_multivariate_and_preserves_target_units(self):
        rng=np.random.default_rng(4);maps=rng.normal(size=(80,5));x=rng.normal(size=(100,5));x[:,1]+=.8*x[:,0]
        operator,tc,scale,condition=e.temporal_operator(e.weights(maps),maps@x.T)
        target=rng.normal(size=(5,12));values=tc@target
        np.testing.assert_allclose(operator@values,target,atol=1e-11)
        np.testing.assert_allclose(tc.std(axis=0,ddof=1),1)
        self.assertGreater(condition,1)

    def test_group_tests_average_runs_before_participants_and_correct_all_60_tests(self):
        rng=np.random.default_rng(19);effects=[];known={}
        for mode in g.inputs.MODES:
            for roi in g.core.ROIS:
                for task in g.core.TASKS:
                    for k in range(1,6):
                        values=rng.normal(.1,.5,132);known[mode,roi,task,k]=values
                        for j,value in enumerate(values):
                            # Unequal run count; participant means must remain the supplied values.
                            for run,delta in ((1,0.),) if j%2 else ((1,-.2),(2,.2)):
                                effects.append(dict(mode=mode,roi=roi,task=task,component=k,subject=str(j),run=run,effect=float(value+delta)))
        participants,groups=e.group_task_tests(effects)
        self.assertEqual(len(groups),60);self.assertEqual(len(participants),60*132)
        for row in groups:
            values=known[row['mode'],row['roi'],row['task'],row['component']]
            self.assertEqual(row['n'],132);self.assertAlmostEqual(row['effect'],values.mean())
            self.assertAlmostEqual(row['p'],stats.ttest_1samp(values,0).pvalue)
            self.assertAlmostEqual(row['p_bonferroni_60'],min(1.,row['p']*60))

    def test_missing_participant_cannot_produce_complete_cohort_task_results(self):
        with self.assertRaisesRegex(ValueError,'Incomplete'):e.group_task_tests([])

    def test_wholebrain_run_matches_direct_residual_regression_for_all_six_bases(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);feat=root/'feat';(feat/'stats').mkdir(parents=True)
            source=root/'source';(source/'run_arrays').mkdir(parents=True);folder=root/'result';folder.mkdir()
            rng=np.random.default_rng(31);shape=(6,6,6);t=60
            reference=nib.Nifti1Image(np.full(shape,1000,dtype='float32'),np.eye(4))
            left=np.zeros(shape,bool);left[:2]=True;right=np.zeros_like(left);right[4:]=True
            masks=dict(left=left,right=right,bilateral=left|right);union=masks['bilateral'];wb=np.ones(shape,bool)
            design=rng.normal(size=(t,4));design-=design.mean(axis=0);vest(feat/'design.mat',design)
            design=g.core.fsl_matrix(feat/'design.mat')
            beta=rng.normal(size=(*shape,4)).astype('float32');native=rng.normal(size=(*shape,t)).astype('float32')
            data=(1000+np.einsum('...p,tp->...t',beta,design)+rng.normal(size=(*shape,t))).astype('float32')
            for name,values in [('filtered_func_data.nii.gz',data),('stats/res4d.nii.gz',native)]:
                image=nib.Nifti1Image(values,reference.affine);image.header.set_xyzt_units('mm','sec');image.header.set_zooms((1,1,1,1.710));nib.save(image,feat/name)
            cope=rng.normal(size=shape).astype('float32');nib.save(nib.Nifti1Image(cope,reference.affine),feat/'stats/cope11.nii.gz')
            for k in range(4):nib.save(nib.Nifti1Image(beta[...,k],reference.affine),feat/f'stats/pe{k+1}.nii.gz')
            joint=data.astype(float)-data.mean(axis=3,keepdims=True,dtype=float)-np.einsum('...p,tp->...t',beta.astype(float),design)
            row=dict(key='test',active_feat=str(feat),filtered_data=str(feat/'filtered_func_data.nii.gz'),native_residual=str(feat/'stats/res4d.nii.gz'),
                subject='10817',task='ugr',run=1,contrast_index=11,nvol=t)
            np.save(source/'run_arrays/test_joint_bold.npy',joint[union].astype('float32'));np.save(source/'run_arrays/test_native_film.npy',native[union])
            primary={}
            for mode in g.inputs.MODES:
                for roi in g.core.ROIS:
                    maps=rng.normal(size=(int(masks[roi].sum()),5));primary[mode,roi]=dict(mask=masks[roi],weights=e.weights(maps))
            result=e.analyze_run(row,source,folder,primary,union,wb,reference)
            self.assertEqual(len(result['effects']),30)
            for mode in g.inputs.MODES:
                residual=joint if mode=='joint_bold' else native
                insula=np.load(source/'run_arrays'/f'test_{mode}.npy')
                for roi in g.core.ROIS:
                    operator,_,_,_=e.temporal_operator(primary[mode,roi]['weights'],insula[masks[roi][union]])
                    expected=(operator@residual[wb].T).T
                    actual=nib.load(folder/f'{mode}_{roi}_connectivity.nii.gz').get_fdata()[wb]
                    np.testing.assert_allclose(actual,expected,atol=1e-6,rtol=1e-5)

    def test_snapshot_requires_all_six_primary_maps_and_preserves_the_running_state(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);analysis=root/'derivatives/insula_ica/unsmoothed-group-ica-fixture';analysis.mkdir(parents=True)
            (analysis/'masks').mkdir();source=root/'record';source.mkdir();fits={}
            for roi in g.core.ROIS:
                mask=analysis/'masks'/(roi+'.nii.gz');mask.write_bytes(b'fixed mask')
                for mode in g.inputs.MODES:
                    folder=analysis/mode/roi;folder.mkdir(parents=True)
                    maps=folder/'oriented_IC.nii.gz';maps.write_bytes(b'completed component map')
                    item=dict(mode=mode,roi=roi,dimensionality=5,state='complete',oriented_maps=str(maps),argv=['melodic','frozen'])
                    marker=dict(item,map_identity=g.core.identity(maps),mask_sha256=g.core.sha(mask),variance_normalization=True)
                    (folder/'complete.json').write_text(json.dumps(marker));fits[f'{mode}/{roi}/dim-5']=item
            rows=[dict(subject=str(j%132),key=str(j)) for j in range(490)]
            (source/'accepted-inputs.json').write_text(json.dumps(rows))
            state=dict(state='running_6_serial_automatic_ICA_fits',output_root=str(analysis),fits=fits)
            path=source/'state.json';path.write_text(json.dumps(state));before=path.read_bytes()
            with patch.object(e,'ROOT',root):
                _,loaded,complete=e.snapshot(source);self.assertEqual(len(complete),6)
                self.assertEqual(path.read_bytes(),before)
                del state['fits']['native_film/right/dim-5'];path.write_text(json.dumps(state))
                with self.assertRaisesRegex(ValueError,'all six|All six'):e.snapshot(source)


if __name__=='__main__':unittest.main()
