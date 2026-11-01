import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock,patch

import nibabel as nib
import numpy as np
import run as f


class FunctionalTests(unittest.TestCase):
    def test_raw_cache_is_reused_and_source_or_content_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);values=np.random.default_rng(4).normal(size=(127,50)).astype('float32')+10000
            loader=Mock(return_value=values)
            actual,reused=f.cache_values(folder,{'source':'first'},loader);self.assertFalse(reused)
            actual,reused=f.cache_values(folder,{'source':'first'},loader);self.assertTrue(reused);self.assertEqual(loader.call_count,1)
            np.testing.assert_array_equal(actual,values)
            with self.assertRaisesRegex(ValueError,'changed'):f.cache_values(folder,{'source':'changed'},loader)
            f.atomic_array(folder/'task_BOLD.npy',values+1)
            with self.assertRaisesRegex(ValueError,'changed'):f.cache_values(folder,{'source':'first'},loader)

    def test_spatial_regression_recovers_sources_including_shared_mean(self):
        rng=np.random.default_rng(32)
        for dim in (3,5,7,11):
            maps=rng.normal(size=(125,dim));sources=rng.normal(size=(dim,90));offset=14000+rng.normal(size=90)*200
            np.testing.assert_allclose(f.spatial_weights(maps)@(maps@sources+offset),sources,atol=2e-10)
        with self.assertRaisesRegex(ValueError,'Nonestimable'):f.spatial_weights(np.ones((127,3)))

    def test_task_fit_uses_standardized_full_design_and_reuses_complete_fit(self):
        rng=np.random.default_rng(11);x=rng.normal(size=(90,8));x[:,6]=1;x[:,7]=0
        responses=(2*x[:,0]-.7*x[:,3]+rng.normal(size=90)*.1)[:,None]+14000
        cs=np.zeros((2,8));cs[0,0]=1;cs[1,3]=1
        row=dict(subject='sample',task='ugr',run=1)
        layout=[dict(mode='functional_mean',roi=f.ROI,dimensionality=0,series='functional_mean',component=0)]
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);effects,reused=f.fit_response(row,folder,{'s':1},responses,layout,x,cs,['a','b']);self.assertFalse(reused)
            zx,_,_,_=f.n.zscore_columns(x,True);zy,_,_,_=f.n.zscore_columns(responses)
            expected=cs@np.linalg.pinv(zx,rcond=1e-10)@zy
            np.testing.assert_allclose([r['effect'] for r in effects],expected[:,0])
            with patch.object(f.n,'fit_timecourses',side_effect=AssertionError('Should not refit')):
                other,reused=f.fit_response(row,folder,{'s':1},responses,layout,x,cs,['a','b'])
            self.assertTrue(reused);self.assertEqual(other,effects)
            with self.assertRaisesRegex(ValueError,'changed'):f.fit_response(row,folder,{'s':2},responses,layout,x,cs,['a','b'])

    def test_written_4d_blocks_are_physically_masked_and_keep_corrected_TR(self):
        rng=np.random.default_rng(5);shape=(6,6,4);reference=nib.Nifti1Image(np.ones(shape),np.eye(4))
        original=np.zeros(np.prod(shape),bool);original[:127]=True;original=original.reshape(shape)
        mask=original.copy();mask.flat[[1,9]]=False
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'masks').mkdir();(root/'run_arrays').mkdir();rows=[]
            f.g.freeze_nifti(root/'masks'/(f.ROI+'.nii.gz'),mask.astype('uint8'),reference.affine)
            for subject in range(132):
                for task in ('trust','ugr'):
                    key=f'sub-{subject}_task-{task}';tr=1.710 if (subject,task)==(0,'ugr') else 1.615
                    row=dict(key=key,subject=str(subject),task=task,run=1,tr=tr);rows.append(row)
                    (root/'run_arrays'/(key+'_complete.json')).write_text('{}')
                    for mode in f.g.inputs.MODES:np.save(root/'run_arrays'/(key+'_'+mode+'.npy'),rng.normal(size=(127,10)).astype('float32'))
            slices=f.build_sessions(rows,root,mask,original,reference)
            for mode in f.g.inputs.MODES:
                folder=root/'sessions'/mode/f.ROI;self.assertEqual(len((folder/'inputs.txt').read_text().splitlines()),264)
                path=folder/'sub-0_task-ugr.nii.gz';im=nib.load(path);data=im.get_fdata()
                self.assertAlmostEqual(im.header.get_zooms()[3],1.710,places=5)
                self.assertEqual(np.count_nonzero(data[~mask[slices]]),0)
                self.assertTrue(np.all(data[mask[slices]].std(axis=1)>0))
            def reject(*args,**kwargs):raise AssertionError('Completed balanced inputs must not be rewritten')
            with patch.object(f.g,'write_nifti',side_effect=reject):f.build_sessions(rows,root,mask,original,reference)
            argv=f.g.melodic_command('joint_bold',f.ROI,0,root,variance_normalization=True,preserve_timing=True)
            self.assertEqual(argv[argv.index('-a')+1],'concat');self.assertNotIn('--vn',argv)
            self.assertNotIn('--tr=1.615',argv);self.assertIn('--disableMigp',argv)

    def test_retained_fit_residuals_use_existing_PEs_and_never_recover_FILM(self):
        rng=np.random.default_rng(6);shape=(6,6,4);t=30
        reference=nib.Nifti1Image(np.ones(shape,dtype='float32'),np.eye(4))
        original=np.zeros(np.prod(shape),bool);original[:127]=True;original=original.reshape(shape)
        design=rng.normal(size=(t,4));y=(14000+rng.normal(size=(*shape,t))*20).astype('float32');pe=rng.normal(size=(*shape,4)).astype('float32')
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);feat=base/'feat';(feat/'stats').mkdir(parents=True);root=base/'out';root.mkdir();(root/'run_arrays').mkdir();record=base/'record';(record/'runs').mkdir(parents=True);parent=base/'parent';(parent/'runs').mkdir(parents=True)
            def write(path,data):nib.save(nib.Nifti1Image(np.asarray(data,dtype='float32'),np.eye(4)),path)
            write(feat/'filtered.nii.gz',y);write(feat/'stats/res4d.nii.gz',rng.normal(size=y.shape));write(feat/'mask.nii.gz',np.ones(shape));write(feat/'stats/sigmasquareds.nii.gz',np.ones(shape))
            for k in range(4):write(feat/f'stats/pe{k+1}.nii.gz',pe[...,k])
            (feat/'design.fsf').write_text('fixture');before=[f.g.core.identity(p,False) for p in feat.rglob('*.nii.gz')]
            row=dict(key='fixture',active_feat=str(feat),filtered_data=str(feat/'filtered.nii.gz'),native_residual=str(feat/'stats/res4d.nii.gz'))
            (parent/'runs/fixture.json').write_text(json.dumps(dict(source_identities=before)))
            with patch.object(f.g.inputs,'inspect_source',return_value=(design,before)),patch.object(f.g.inputs,'ensure_native',side_effect=AssertionError('FILM recovery forbidden')):
                qc,reused=f.prepare_residuals(row,root,record,original,reference,parent)
                self.assertFalse(reused);self.assertFalse(qc['native_recovered_by_film'])
                expected=y[original].astype(float);expected-=expected.mean(axis=1,keepdims=True);expected-=pe[original].astype(float)@design.T
                np.testing.assert_allclose(np.load(root/'run_arrays/fixture_joint_bold.npy'),expected,rtol=1e-6,atol=1e-4)
                _,reused=f.prepare_residuals(row,root,record,original,reference,parent);self.assertTrue(reused)
            changed=[dict(v,size=v['size']+1) for v in before]
            with patch.object(f.g.inputs,'inspect_source',return_value=(design,changed)):
                with self.assertRaisesRegex(ValueError,'differ'):f.prepare_residuals(row,root,record,original,reference,parent)

    def test_existing_all_IC_barplots_get_five_bar_mean_without_changing_source_rows(self):
        source=f.ROOT/'code/insula_ica/records/linux2-task-dimensions-20261030T034620Z-0c19a6/group_task_effects.tsv'
        rows=f.g.core.read_table(source);before=json.dumps(rows,sort_keys=True)
        means=[dict(series='functional_mean',contrast=label,effect=.1,ci_low=.05,ci_high=.15) for label in f.n.BAR_LABELS]
        with tempfile.TemporaryDirectory() as tmp:
            f.insula_plots(rows,means,Path(tmp));self.assertEqual(len(list(Path(tmp).glob('*.png'))),7)
        self.assertEqual(json.dumps(rows,sort_keys=True),before)


if __name__=='__main__':unittest.main()
