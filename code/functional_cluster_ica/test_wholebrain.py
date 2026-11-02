import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np
import wholebrain as w
import network_summary as s


def image(path,data,tr=None):
    im=nib.Nifti1Image(np.asarray(data,dtype='float32'),np.eye(4))
    if tr is not None:im.header.set_zooms((1,1,1,tr))
    nib.save(im,path)


class WholebrainTests(unittest.TestCase):
    def test_wholebrain_stage2_matches_direct_GLS_residual_regression_at_both_dimensions(self):
        rng=np.random.default_rng(51);shape=(6,6,6);t=60;mask=np.zeros(np.prod(shape),bool);mask[:117]=1;mask=mask.reshape(shape);original=np.zeros(np.prod(shape),bool);original[:127]=1;original=original.reshape(shape);wb=np.ones(shape,bool)
        design=rng.normal(size=(t,4));design-=design.mean(axis=0);pe=rng.normal(size=(*shape,4)).astype('float32');native=rng.normal(size=(*shape,t)).astype('float32')
        data=(12000+rng.normal(size=(*shape,t))*10).astype('float32');joint=(data.astype(float)-data.mean(axis=3,keepdims=True,dtype=float))-np.einsum('...p,tp->...t',pe.astype(float),design)
        reference=nib.Nifti1Image(np.ones(shape),np.eye(4))
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);feat=root/'feat';(feat/'stats').mkdir(parents=True);source=root/'source';(source/'run_arrays').mkdir(parents=True);out=root/'out';out.mkdir()
            w.vest(feat/'design.mat',design);design=w.f.g.core.fsl_matrix(feat/'design.mat')
            image(feat/'filtered_func_data.nii.gz',data,1.710);image(feat/'stats/res4d.nii.gz',native,1.710)
            for k in range(4):image(feat/f'stats/pe{k+1}.nii.gz',pe[...,k])
            # Recompute with the actual written design precision.
            joint=(data.astype(float)-data.mean(axis=3,keepdims=True,dtype=float))-np.einsum('...p,tp->...t',pe.astype(float),design)
            for mode,values in [('joint_bold',joint),('native_film',native)]:np.save(source/'run_arrays'/f'fixture_{mode}.npy',values[original].astype('float32'))
            families={}
            for mode in w.f.g.inputs.MODES:
                for dim in (3,5):
                    key=f'{mode}_{dim}';maps=rng.normal(size=(117,dim));families[key]=dict(mode=mode,dim=dim,weights=w.f.spatial_weights(maps))
            row=dict(key='fixture',active_feat=str(feat),filtered_data=str(feat/'filtered_func_data.nii.gz'),native_residual=str(feat/'stats/res4d.nii.gz'),nvol=t)
            w.analyze_run(row,source,out,families,original,mask,wb,reference)
            for key,item in families.items():
                residual=np.load(source/'run_arrays'/f'fixture_{item["mode"]}.npy')[mask[original]];operator,tc,_=w.temporal_operator(item['weights'],residual)
                target=joint if item['mode']=='joint_bold' else native.astype(float);target=target-target.mean(axis=3,keepdims=True)
                expected=(operator@target[wb].T).T;actual=nib.load(out/(key+'_DR.nii.gz')).get_fdata()[wb]
                np.testing.assert_allclose(actual,expected,rtol=1e-5,atol=2e-6)
                np.testing.assert_allclose(tc.mean(axis=0),0,atol=1e-12);np.testing.assert_allclose(tc.std(axis=0,ddof=1),1)

    def test_task_collapses_are_paired_inside_participant_before_CI(self):
        rows=[]
        for k in range(132):
            base=float(k);v=(base,base+.2,base,base+.7,base+.1)
            for label,value in zip(w.f.n.BAR_LABELS,v):rows.append(dict(family='fixture',subject=str(k),contrast=label,effect=value))
        participants,groups=s.paired_summaries(rows,('family',));ug=next(r for r in groups if r['contrast']=='UGR social − nonsocial');tg=next(r for r in groups if r['contrast']=='Trust stranger − computer')
        self.assertAlmostEqual(ug['effect'],.2);self.assertAlmostEqual(tg['effect'],.1);self.assertLess(ug['se'],1e-13);self.assertLess(tg['se'],1e-13)
        with self.assertRaisesRegex(ValueError,'Duplicate'):s.paired_summaries(rows+[rows[0]],('family',))
        with self.assertRaisesRegex(ValueError,'Full paired'):s.paired_summaries(rows[:-1],('family',))

    def test_participant_maps_equal_weight_runs_then_tasks_and_cache_merged_inputs(self):
        reference=nib.Nifti1Image(np.ones((2,2,2)),np.eye(4));wb=np.ones((2,2,2),bool);key='fixture';item=dict(key=key,arm='sm0',dim=3)
        rows=[];markers={}
        values=[('a','ugr',1,2),('a','ugr',2,4),('a','trust',1,9),('b','ugr',1,1),('b','trust',1,2),('b','trust',2,6)]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);record=root/'record';(record/'group_maps').mkdir(parents=True)
            for sub,task,run,value in values:
                run_key=f'{sub}_{task}_{run}';rows.append(dict(key=run_key,subject=sub,task=task,run=run));markers['sm0',run_key]=run_key
                folder=root/'runs/sm0'/run_key;folder.mkdir(parents=True);image(folder/(key+'_DR.nii.gz'),np.full((2,2,2,3),value))
            subjects=w.aggregate_family(item,{'sm0':rows},root,record,reference,wb,markers);self.assertEqual(subjects,['a','b'])
            pooled=nib.load(record/'group_maps/fixture_pooled_mean.nii.gz').get_fdata();np.testing.assert_allclose(pooled,4.25)
            merged=root/'merged/fixture/IC1.nii.gz';im=nib.load(merged).get_fdata();np.testing.assert_allclose(im[...,0],6);np.testing.assert_allclose(im[...,1],2.5)
            before=w.f.g.core.identity(merged);w.aggregate_family(item,{'sm0':rows},root,record,reference,wb,markers);self.assertEqual(w.f.g.core.identity(merged),before)

    def test_group_design_is_frozen_and_randomise_has_two_explicit_signed_contrasts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);d=root/'group.mat';c=root/'group.con';w.vest(d,np.ones((132,1)));w.vest(c,[[1],[-1]],True)
            before=w.f.g.core.identity(d);w.vest(d,np.ones((132,1)));self.assertEqual(w.f.g.core.identity(d),before)
            np.testing.assert_array_equal(w.f.g.core.fsl_matrix(c),[[1],[-1]])
            argv=w.randomise_command(root/'merged.nii.gz',root/'out',root/'mask.nii.gz',d,c,5000)
            self.assertIn('-T',argv);self.assertNotIn('-D',argv);self.assertEqual(argv[argv.index('-n')+1],'5000')
            with self.assertRaisesRegex(ValueError,'changed'):w.vest(d,np.ones((132,1))*2)

    def test_completed_permutation_checkpoint_is_reused_without_external_call(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);folder=root/'result';folder.mkdir();merged=root/'merged.nii.gz';mask=root/'mask.nii.gz';d=root/'group.mat';c=root/'group.con'
            image(merged,np.ones((2,2,2,8)));image(mask,np.ones((2,2,2)));w.vest(d,np.ones((8,1)));w.vest(c,[[1],[-1]],True)
            command=w.randomise_command(merged,folder/'dr',mask,d,c,5000)
            sig=dict(command=command,merged=w.f.g.core.identity(merged),mask=w.f.g.core.identity(mask),design=w.f.g.core.identity(d),contrasts=w.f.g.core.identity(c),family_tests=64,components=3)
            output=folder/'result.nii.gz';image(output,np.ones((2,2,2)));w.f.save(folder/'complete.json',dict(signature=sig,outputs=[w.f.g.core.identity(output)]))
            with patch.object(w.subprocess,'run',side_effect=AssertionError('Completed permutations must not rerun')):r=w.randomise_job(merged,folder,mask,d,c,5000,3)
            self.assertTrue(r['reused'])


if __name__=='__main__':unittest.main()
