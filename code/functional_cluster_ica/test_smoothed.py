import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np
import smoothed as m


class MatchedTests(unittest.TestCase):
    def fixture(self,root):
        shape=(4,4,3);t=24;ref=nib.Nifti1Image(np.ones(shape,dtype='float32'),np.eye(4))
        x=np.random.default_rng(7).normal(size=(t,4));nuisance=x[:,-1:]
        bold=root/'raw.nii.gz';im=nib.Nifti1Image(np.ones((*shape,t),dtype='float32'),np.eye(4));im.header.set_zooms((1,1,1,1.710));nib.save(im,bold)
        conf=root/'nuisance.txt';np.savetxt(conf,nuisance);ev=root/'task.txt';np.savetxt(ev,[[1,0,1],[3,0,1]])
        paths=[]
        for name,smooth in [('sm0',0),('sm5',5)]:
            feat=root/name;feat.mkdir();(feat/'stats').mkdir();paths.append(feat)
            settings={'fmri(smooth)':smooth,'fmri(outputdir)':str(feat),'fmri(overwrite_yn)':0,'fmri(tr)':1.710,'fmri(custom1)':str(ev),'fmri(shape1)':3,'fmri(temphp_yn)':0,'feat_files(1)':str(bold),'confoundev_files(1)':str(conf)}
            (feat/'design.fsf').write_text('\n'.join(f'set {k} "{v}"' for k,v in settings.items())+'\n')
            for filename,values in [('design.mat',x),('design.con',np.eye(4))]:
                with (feat/filename).open('w') as out:out.write('/Matrix\n');np.savetxt(out,values)
            np.savetxt(feat/'confoundevs.txt',nuisance);nib.save(im,feat/'stats/res4d.nii.gz')
        row=dict(key='sub-10817_ses-01_task-ugr_run-1',active_feat=str(paths[0]),confounds=str(conf),nvol=t,tr=1.710,task='ugr',confounds_sha256=m.f.g.core.sha(conf))
        signature=dict(scientific_settings=m.scientific_settings(m.f.g.core.fsl_settings(paths[0]/'design.fsf')),BOLD=m.f.g.core.identity(bold,False),TR=1.710)
        return ref,row,signature,paths

    def test_only_smoothing_and_output_paths_change_in_generated_FSF(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref,row,signature,paths=self.fixture(root);target=root/'new.feat'
            before=m.f.g.core.fsl_settings(paths[0]/'design.fsf')
            text=m.smoothed_fsf((paths[0]/'design.fsf').read_text(),row,target);p=root/'new.fsf';p.write_text(text);after=m.f.g.core.fsl_settings(p)
            self.assertEqual(after['fmri(smooth)'],'5');self.assertEqual(after['fmri(outputdir)'],str(target))
            self.assertEqual(m.scientific_settings(before),m.scientific_settings(after))
            self.assertEqual(after['feat_files(1)'],before['feat_files(1)']);self.assertEqual(after['fmri(tr)'],'1.71')

    def test_matching_fit_accepts_same_BOLD_design_nuisance_and_rejects_each_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);reference,row,signature,paths=self.fixture(root);feat=paths[1];original=(feat/'design.fsf').read_text()
            with patch.object(m.f.g.inputs,'inspect_source',return_value=(np.zeros((24,4)),[])):
                selected,_=m.check_matched(row,feat,signature,reference);self.assertEqual(selected['smoothing_mm'],5)
                for old,new in [('fmri(tr) "1.71"','fmri(tr) "1.615"'),('fmri(temphp_yn) "0"','fmri(temphp_yn) "1"')]:
                    self.assertIn(old,original);(feat/'design.fsf').write_text(original.replace(old,new))
                    with self.assertRaises(ValueError):m.check_matched(row,feat,signature,reference)
                (feat/'design.fsf').write_text(original)
                other=root/'other.nii.gz';other.write_bytes((root/'raw.nii.gz').read_bytes());(feat/'design.fsf').write_text(original.replace(str(root/'raw.nii.gz'),str(other)))
                with self.assertRaisesRegex(ValueError,'different BOLD'):m.check_matched(row,feat,signature,reference)
                (feat/'design.fsf').write_text(original);p=feat/'design.mat';saved=p.read_bytes();p.write_text('/Matrix\n'+'\n'.join('0 0 0 0' for _ in range(24)))
                with self.assertRaisesRegex(ValueError,'temporal model'):m.check_matched(row,feat,signature,reference)
                p.write_bytes(saved);np.savetxt(feat/'confoundevs.txt',np.zeros((24,1)))
                with self.assertRaisesRegex(ValueError,'nuisance'):m.check_matched(row,feat,signature,reference)

    def test_completed_FEAT_checkpoint_prevents_any_FSL_refit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref,row,signature,paths=self.fixture(root);record=root/'record';(record/'FEAT-snapshot/runs').mkdir(parents=True);(root/'FEAT-checkpoints').mkdir()
            ids=[m.f.g.core.identity(paths[1]/'design.mat')]
            m.f.save(root/'FEAT-checkpoints'/(row['key']+'.json'),dict(signature=signature,selected_feat=str(paths[1]),outputs=ids,origin='verified_existing_sm5'))
            with patch.object(m,'source_signature',return_value=(signature,{},None)),patch.object(m,'check_matched',return_value=(dict(row,active_feat=str(paths[1])),ids)),patch.object(m.subprocess,'run',side_effect=AssertionError('No FEAT or feat_model should run')):
                selected,status=m.fit_feat(row,root,record,root,ref,[])
            self.assertTrue(status['reused_completed_feat']);self.assertEqual(selected['active_feat'],str(paths[1]))

    def test_wrong_regenerated_model_stops_before_FEAT_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);reference,row,signature,paths=self.fixture(base);root=base/'new';root.mkdir();record=base/'record';record.mkdir()
            for name in ('fsf','logs','feat','FEAT-checkpoints'):(root/name).mkdir()
            calls=[]
            def model_only(argv,**kwargs):
                calls.append(argv[0])
                self.assertEqual(argv[0],'feat_model','A mismatched design must stop before FEAT')
                stem=Path(argv[1])
                for suffix,data in (('.mat',np.zeros((24,4))),('.con',np.eye(4))):
                    with stem.with_suffix(suffix).open('w') as stream:stream.write('/Matrix\n');np.savetxt(stream,data)
            with patch.object(m,'source_signature',return_value=(signature,{},None)),patch.object(m.subprocess,'run',side_effect=model_only):
                with self.assertRaisesRegex(ValueError,'before FEAT'):m.fit_feat(row,root,record,record,reference,[])
            self.assertEqual(calls,['feat_model'])

    def test_raw_task_cache_is_reused_without_any_image_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);key='fixture';folder=root/'task_values'/key;values=np.random.default_rng(4).normal(size=(127,30)).astype('float32');anchors=dict(raw='unchanged')
            m.f.cache_values(folder,dict(anchors,mask_sha256=m.f.MASK_SHA),lambda:values)
            with patch.object(m.f,'task_context',return_value=(np.eye(30),np.eye(30),['fixture'],None,anchors)),patch.object(m.f,'verify_anchors'),patch.object(m.nib,'load',side_effect=AssertionError('No BOLD image should be read')):
                actual,*_=m.load_task_cache(dict(key=key),root,None,root)
            np.testing.assert_array_equal(actual,values)

    def test_cross_smoothing_comparisons_match_maps_before_paired_task_differences(self):
        rng=np.random.default_rng(9);shape=(4,4,2);mask=np.ones(shape,bool);a=rng.normal(size=(*shape,3));b=np.stack([-a[...,2],a[...,0],a[...,1]],axis=3)
        key='joint_bold/functional_cluster/dim-3';prior=[];after=[]
        for subject in range(132):
            for contrast in m.f.n.BAR_LABELS:
                values=rng.normal(size=3)
                for k,value in enumerate(values,1):prior.append(dict(subject=str(subject),contrast=contrast,mode='joint_bold',dimensionality=3,component=k,effect=float(value)))
                for k,value in enumerate((-(values[2]+.15),values[0]+.15,values[1]+.15),1):after.append(dict(subject=str(subject),contrast=contrast,mode='joint_bold',dimensionality=3,component=k,effect=float(value)))
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';source.mkdir();record=root/'record';record.mkdir()
            for name,data in [('sm0',a),('sm5',b)]:nib.save(nib.Nifti1Image(data,np.eye(4)),root/(name+'.nii.gz'))
            first=dict(state='complete',oriented_maps=str(root/'sm0.nii.gz'));second=dict(state='complete',oriented_maps=str(root/'sm5.nii.gz'),mode='joint_bold',dimensionality=3)
            m.f.save(source/'state.json',dict(fits={key:first}));m.f.g.core.write_table(source/'functional_IC_participant_effects.tsv',prior)
            m.compare_smoothing(source,dict(fits={key:second}),mask,record,after)
            result=m.f.g.core.read_table(record/'smoothing_paired_task_comparisons.tsv');self.assertEqual(len(result),15)
            np.testing.assert_allclose([float(r['paired_sm5_minus_sm0']) for r in result],.15,atol=1e-12)
            np.testing.assert_allclose([abs(float(r['spatial_r'])) for r in result],1,atol=1e-12)


if __name__=='__main__':unittest.main()
