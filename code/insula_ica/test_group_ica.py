from pathlib import Path
from contextlib import ExitStack
import json
import os
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

import group_ica_inputs as i
import run_group_ica as g


class PairedGroupTests(unittest.TestCase):
    def test_bold_residual_uses_supplied_gls_coefficients_without_refit(self):
        rng=np.random.default_rng(12);t=80
        design=np.column_stack([np.sin(np.arange(t)/8),np.cos(np.arange(t)/13)])
        design-=design.mean(axis=0)
        data=rng.normal(size=(7,t))+500
        beta=rng.normal(size=(7,2))*4
        residual=i.joint_bold_residual(data,design,beta)
        np.testing.assert_allclose(residual+beta@design.T,data-data.mean(axis=1,keepdims=True))
        ols=(data-data.mean(axis=1,keepdims=True))@np.linalg.pinv(design).T
        self.assertGreater(np.linalg.norm(beta-ols),1)
        np.testing.assert_allclose(residual,i.joint_bold_residual(data+123,design,beta),atol=1e-10)

    def test_decomposition_rejects_missing_or_nonfinite_coefficients(self):
        with self.assertRaises(ValueError):i.joint_bold_residual(np.ones((3,5)),np.ones((5,2)),np.ones((3,1)))
        with self.assertRaises(ValueError):i.joint_bold_residual(np.full((3,5),np.nan),np.ones((5,2)),np.ones((3,2)))

    def test_real_record_selection_retains_flagged_runs_and_repaired_pair(self):
        records=Path(__file__).parent/'records';rows=i.select_trust(records)
        self.assertEqual(len(rows),241)
        bykey={r['key']:r for r in rows}
        for key in i.ACCEPTED_FLAGS:
            self.assertIn(i.STOCK,bykey[key]['source_selection_record'])
        for run in (1,2):
            row=bykey[f'sub-10834_ses-01_task-trust_run-{run}']
            self.assertIn(i.REPAIR,row['source_selection_record'])
            self.assertIn('registration-sub-10834-ignore-sbref',row['native_residual'])
        self.assertIn(i.BASE,bykey['sub-10369_ses-01_task-trust_run-1']['source_selection_record'])

    def test_full_cohort_selection_uses_frozen_ugr_flags_not_directory_globs(self):
        records=Path(__file__).parent/'records'
        trust=i.select_trust(records)
        cohort=sorted({r['subject'] for r in trust})
        ugr=[]
        for j,subject in enumerate(cohort):
            for run in ((1,2) if j<117 else (1,)):
                ugr.append(dict(subject=subject,task='ugr',status='included',run=str(run),nvol='240',tr='1.615',
                    key=f'sub-{subject}_ses-01_task-ugr_run-{run}',original_feat=f'/original/sub-{subject}/run-{run}.feat',
                    filtered_data=f'/original/sub-{subject}/run-{run}.feat/filtered_func_data.nii.gz',confound_matrix='/confounds'))
        ugr.append(dict(task='ugr',status='excluded_original_qc',subject='outside'))
        with patch.object(i.core,'read_table',return_value=ugr),patch.object(i.core,'sha',return_value='frozen'):
            rows=i.select_inputs(Path('/frozen-audit'))
        self.assertEqual(len(rows),490)
        self.assertEqual(sum(r['task']=='ugr' for r in rows),249)
        self.assertNotIn('outside',{r['subject'] for r in rows})

    def test_balancing_preserves_all_frames_and_equal_run_energy(self):
        rng=np.random.default_rng(2)
        a=rng.normal(size=(9,40))*8;b=rng.normal(size=(9,70))*.2
        values,weights=i.balanced_session([a,b])
        self.assertEqual(values.shape,(9,110))
        np.testing.assert_allclose(np.sum(values[:,:40]**2,axis=1),.5,rtol=1e-6)
        np.testing.assert_allclose(np.sum(values[:,40:]**2,axis=1),.5,rtol=1e-6)
        self.assertEqual(len(weights),2)

    def test_crop_preserves_mni_coordinates_exactly(self):
        affine=np.diag([2.7,2.7,2.97,1]);affine[:3,3]=[-75,-110,-72]
        reference=nib.Nifti1Image(np.zeros((20,20,20)),affine)
        mask=np.zeros(reference.shape,dtype=bool);mask[3:9,4:11,5:12]=True
        slices,cropped=i.crop_geometry(reference,mask)
        origin=np.array([s.start for s in slices]);point=np.array([2,3,4])
        np.testing.assert_allclose(nib.affines.apply_affine(cropped,point),nib.affines.apply_affine(affine,origin+point))

    def test_plan_has_all_24_fits_with_fixed_and_automatic_dimensions(self):
        plans=g.plans();self.assertEqual(len(plans),len(set(plans)))
        self.assertEqual(len(plans),24)
        for mode in i.MODES:
            for roi in g.core.ROIS:self.assertEqual({d for m,r,d in plans if m==mode and r==roi},{0,3,5,7})
        argv=g.melodic_command('joint_bold','left',5,Path('/analysis'))
        self.assertIn('--seed=20261028',argv)
        self.assertIn('--vn',argv);self.assertIn('--update_mask',argv)
        self.assertEqual(argv[argv.index('-i')+1],'/analysis/sessions/joint_bold/left/inputs.txt')
        self.assertNotIn('--filter',argv);self.assertNotIn('--bgthreshold',argv)

    def test_matching_recovers_permutation_and_signs(self):
        rng=np.random.default_rng(9);a=rng.normal(size=(80,5))
        permutation=[3,0,4,1,2];b=a[:,permutation]*np.array([-1,1,-1,1,1])
        result=g.compare_maps(a,b)
        self.assertTrue(all(abs(x['absolute_spatial_r']-1)<1e-10 for x in result['matches']))
        self.assertLess(max(result['centered_subspace_angles_degrees']),1e-10)
        self.assertEqual(len(result['all_pair_spatial_r']),5)

    def test_automatic_dimension_comparison_allows_different_component_counts(self):
        a=np.random.default_rng(18).normal(size=(80,5));b=a[:,:3]*[-1,1,1]
        result=g.compare_maps(a,b)
        self.assertEqual(result['joint_bold_components'],5);self.assertEqual(result['native_film_components'],3)
        self.assertEqual(len(result['matches']),3)
        self.assertEqual(result['unmatched_joint_bold_components'],[4,5])
        self.assertTrue(all(abs(x['absolute_spatial_r']-1)<1e-10 for x in result['matches']))

    def test_missing_trust_residual_never_triggers_preprocessing_or_ugr_replay(self):
        row=dict(task='trust',native_residual='/missing/res4d.nii.gz')
        with patch.object(i.subprocess,'run') as run:
            with self.assertRaises(ValueError):i.ensure_native(row,Path('/analysis'))
            run.assert_not_called()

    def test_common_support_is_shared_between_hemispheres_and_both_scales(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);record=root/'record';record.mkdir();(root/'run_arrays').mkdir()
            shape=(4,4,4);left=np.zeros(shape,bool);right=left.copy();left[:2]=True;right[2:]=True
            anatomical=dict(left=left,right=right,bilateral=left|right)
            a=np.ones(64,bool);b=a.copy();b[0]=False
            np.save(root/'run_arrays/a_valid.npy',a);np.save(root/'run_arrays/b_valid.npy',b)
            masks=g.fixed_masks([dict(key='a'),dict(key='b')],root,record,anatomical,nib.Nifti1Image(np.zeros(shape),np.eye(4)))
            self.assertEqual(masks['bilateral'].sum(),63)
            self.assertTrue(np.array_equal(masks['left']|masks['right'],masks['bilateral']))
            self.assertFalse(masks['left'][0,0,0])

    def test_frozen_ugr_inputs_reject_model_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);old=root/'original';old.mkdir();audit=root/'audit';audit.mkdir()
            names=['design.mat','design.con','mask.nii.gz','mean_func.nii.gz','filtered_func_data.nii.gz','confounds.txt']
            for name in names:(old/name).write_bytes(b'frozen historical inputs')
            frozen=[i.core.identity(old/name,content=not name.endswith('.nii.gz')) for name in names]
            (audit/'audit.json').write_text(json.dumps(dict(inputs=frozen)))
            row=dict(task='ugr',original_feat=str(old),filtered_data=str(old/'filtered_func_data.nii.gz'),confounds=str(old/'confounds.txt'))
            i.verify_frozen_ugr([row],audit)
            (old/'design.mat').write_text('changed model')
            with self.assertRaises(ValueError):i.verify_frozen_ugr([row],audit)

    def test_residual_checkpoint_reuses_arrays_without_film_or_reconstruction(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);record=root/'record';record.mkdir();arrays=root/'run_arrays';arrays.mkdir()
            native=root/'native.nii.gz';native.write_bytes(b'unchanged source')
            path=arrays/'run_joint_bold.npy';np.save(path,np.ones((5,20),dtype='float32'))
            qc=dict(source_identities=[],native_residual_identity=i.core.identity(native,False))
            (arrays/'run_complete.json').write_text(json.dumps(dict(qc=qc,arrays=[i.core.identity(path)])))
            with patch.object(i,'inspect_source',return_value=(np.ones((20,2)),[])),patch.object(i,'ensure_native') as replay:
                self.assertEqual(i.prepare_run(dict(key='run'),root,record,None,None),qc)
                replay.assert_not_called()

    def test_end_to_end_orchestration_builds_both_scales_all_24_fits_and_compares(self):
        with tempfile.TemporaryDirectory() as directory,ExitStack() as p:
            base=Path(directory);root=base/'analysis';root.mkdir();record=base/'record';record.mkdir()
            for name in ('figures','errors'):(record/name).mkdir()
            old=base/'retained';old.mkdir();aff=np.eye(4);aff[0,3]=-3
            reference=nib.Nifti1Image(np.full((6,6,6),1000,dtype='float32'),aff)
            nib.save(reference,old/'mean_func.nii.gz');native=old/'res4d.nii.gz';native.touch()
            fsldir=base/'fsl';atlasdir=fsldir/'data/atlases/HarvardOxford';atlasdir.mkdir(parents=True)
            labels=np.full((6,6,6),2,dtype='uint8');labels[3]=0
            nib.save(nib.Nifti1Image(labels,aff),atlasdir/'HarvardOxford-cort-maxprob-thr25-2mm.nii.gz')
            (atlasdir.parent/'HarvardOxford-Cortical.xml').write_text('<atlas><data><label index="1">Insular Cortex</label></data></atlas>')
            rows=[]
            for j in range(132):
                for task in ('trust','ugr'):
                    runs=([1] if j<131 else [])+([2] if j<109 or j==131 else []) if task=='trust' else [1]+([2] if j<117 else [])
                    for run in runs:
                        subject=str(10000+j)
                        rows.append(dict(subject=subject,task=task,run=run,nvol=20,nss=0,
                            key=f'sub-{subject}_ses-01_task-{task}_run-{run}',active_feat=str(old),native_residual=str(native)))
            self.assertEqual(len(rows),490)
            p.enter_context(patch.dict(os.environ,FSLDIR=str(fsldir)))
            p.enter_context(patch.object(g.core,'software',return_value=dict(fsl='test')))
            p.enter_context(patch.object(g.shutil,'disk_usage',return_value=SimpleNamespace(free=1024**4)))
            p.enter_context(patch.object(g.subprocess,'run',return_value=SimpleNamespace(stdout='--update_mask --disableMigp --vn',stderr='')))
            p.enter_context(patch.object(i,'select_inputs',return_value=rows));p.enter_context(patch.object(i,'verify_frozen_ugr'))
            p.enter_context(patch.object(i,'inspect_source'));p.enter_context(patch.object(g,'pilot_checks'))
            def prepare(row,root,record,union,reference):
                folder=root/'run_arrays';folder.mkdir(exist_ok=True);rng=np.random.default_rng(row['run'])
                for mode in i.MODES:np.save(folder/(row['key']+'_'+mode+'.npy'),rng.normal(size=(int(union.sum()),20)).astype('float32'))
                np.save(folder/(row['key']+'_valid.npy'),np.ones(int(union.sum()),bool))
                (folder/(row['key']+'_complete.json')).write_text('{}')
                return dict(native_recovered_by_film=False)
            p.enter_context(patch.object(i,'prepare_run',side_effect=prepare))
            def fit(mode,roi,dim,root,masks,slices,reference):
                folder=root/'ica'/mode/roi/f'dim-{dim}';folder.mkdir(parents=True)
                actual=dim or (4 if mode=='joint_bold' else 6)
                rng=np.random.default_rng(actual);data=np.zeros((*reference.shape,actual),dtype='float32')
                data[masks[roi]]=rng.normal(size=(int(masks[roi].sum()),actual))
                path=folder/'oriented_IC.nii.gz';g.core.save_image(path,data,reference)
                return dict(mode=mode,roi=roi,dimensionality=dim,actual_components=actual,oriented_maps=str(path))
            p.enter_context(patch.object(g,'fit_one',side_effect=fit))
            p.enter_context(patch.object(g,'component_figure',side_effect=lambda maps,mask,bg,affine,path,label,**kwargs:path.write_bytes(b'fixture image')))
            state=dict(runs={},fits={});g.run(SimpleNamespace(audit=base/'audit',jobs=64,ica_jobs=18,auto_ica_jobs=6),root,record,state)
            self.assertEqual(state['state'],'paired_ICA_complete_review_pending')
            self.assertEqual(len(state['runs']),490);self.assertEqual(len(state['fits']),24)
            comparisons=json.loads((record/'between-scale-spatial-comparisons.json').read_text())
            self.assertEqual(len(comparisons),12)
            self.assertEqual(len(list((record/'figures').glob('*.png'))),24)
            for mode in i.MODES:
                for roi in g.core.ROIS:
                    folder=root/'sessions'/mode/roi
                    listing=(folder/'inputs.txt').read_text().splitlines()
                    self.assertEqual(len(listing),264)
                    image=nib.load(listing[0]);data=image.get_fdata()
                    mask=nib.load(root/'masks'/(roi+'_crop.nii.gz')).get_fdata()>0
                    self.assertEqual(np.count_nonzero(data[~mask]),0)
                    self.assertTrue((data[mask].std(axis=1)>0).all())
                    g.verify_ica_inputs(mode,roi,root)
            Path(listing[0]).write_bytes(b'changed after physical-volume verification')
            with self.assertRaisesRegex(ValueError,'Previously verified 4D input changed'):
                g.verify_ica_inputs(mode,roi,root)

    def test_physical_mask_check_rejects_data_in_other_hemisphere(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.nii.gz';mask=np.zeros((4,4,4),bool);mask[:2]=True
            mask_image=nib.Nifti1Image(mask.astype('uint8'),np.eye(4))
            values=np.zeros((4,4,4,20),dtype='float32')
            values[mask]=np.random.default_rng(4).normal(size=(32,20))
            g.write_nifti(path,values,np.eye(4))
            self.assertTrue(g.validate_masked_4d(path,mask_image,20)['physical_4d_mask_verified'])
            values[3,1,1]=1;g.write_nifti(path,values,np.eye(4))
            with self.assertRaisesRegex(ValueError,'outside the specific insula ROI'):
                g.validate_masked_4d(path,mask_image,20)


if __name__=='__main__':unittest.main()
