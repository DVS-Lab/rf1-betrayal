import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

import rerun_activation as r


class ActivationRepairTests(unittest.TestCase):
    def test_input_repair_is_subject_and_contrast_specific(self):
        tr=Path('/corrected/corrected_TR_1p710.feat')
        old='/project/sub-10817/ses-01/L1_task-ugr_model-3_type-act_run-1_sm-5.feat/stats/cope11.nii.gz'
        self.assertEqual(r.choose_input(old,11,tr)[0],tr/'stats/cope11.nii.gz')
        old='/project/sub-11031/ses-01/L2_task-ugr_model-3_type-act_sm-5.gfeat/cope11.feat/stats/cope1.nii.gz'
        self.assertIn('run-2_sm-5.feat/stats/cope11',str(r.choose_input(old,11,tr)[0]))
        unchanged=old.replace('cope11.feat','cope13.feat')
        self.assertEqual(str(r.choose_input(unchanged,13,tr)[0]),unchanged)
        with self.assertRaises(ValueError):r.choose_input(old,12,tr)

    def test_fsf_changes_only_named_settings(self):
        original='set fmri(outputdir) "/old"\nset feat_files(1) "/old/cope"\nset fmri(evg1.1) 0.123\nset fmri(z_thresh) 3.1\n'
        changed=r.replace_settings(original,{'fmri(outputdir)':'/new','feat_files(1)':'/new/cope'})
        self.assertIn('set fmri(evg1.1) 0.123\nset fmri(z_thresh) 3.1',changed)
        with self.assertRaises(ValueError):r.replace_settings(original,{'missing':'x'})

    def test_real_saved_six_group_FSFs_plan_retains_132_inputs_and_all_covariates(self):
        real=r.ROOT; audit=real/'code/ugr_reanalysis/records/linux2-original-group-audit-20261030T113752Z-f372a9'
        original_manifest=json.loads((audit/'manifest.json').read_text())
        copies={v['source']:audit/v['copy'] for v in original_manifest if 'copy' in v}
        cases=[v for v in json.loads((audit/'group_cases.json').read_text()) if v['task']=='ugr']
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); out=root/'output';out.mkdir();rec=root/'record';rec.mkdir()
            for name in ('fsfs','sources'):(rec/name).mkdir()
            ar=root/'audit';ar.mkdir();tr=root/'timing';tr.mkdir();fsldir=root/'fsl';(fsldir/'data/standard').mkdir(parents=True)
            image=nib.Nifti1Image(np.ones((3,4,5),dtype='float32'),np.eye(4))
            nib.save(image,fsldir/'data/standard/MNI152_T1_2mm_brain.nii.gz')
            (root/'code').mkdir();shutil.copyfile(real/'code/sublist_n132.txt',root/'code/sublist_n132.txt')
            prefix='/ZPOOL/data/projects/rf1-betrayal'
            translate=lambda value:str(value).replace(prefix,str(root))
            manifest=[]
            for v in original_manifest:
                if 'copy' not in v:continue
                if '/sub-10817/' not in v['source'] and '/sub-11031/' not in v['source'] and v['role'] not in ('group_design','original_group_map'):continue
                target=Path(translate(v['source']));target.parent.mkdir(parents=True,exist_ok=True)
                if target.suffix=='.fsf':target.write_text(translate(copies[v['source']].read_text()))
                elif target.name.endswith('.nii.gz'):nib.save(image,target)
                else:shutil.copyfile(copies[v['source']],target)
                manifest.append(dict(source=str(target),identity=r.core.identity(target)))
            (ar/'manifest.json').write_text(json.dumps(manifest))
            for name in ('state.json','original_UGR_contrast_audit.tsv'):shutil.copyfile(audit/name,ar/name)
            translated=json.loads(translate(json.dumps(cases)));(ar/'group_cases.json').write_text(json.dumps(translated))
            for case in translated:
                f=Path(case['path'])/'cope1.feat';f.mkdir(parents=True,exist_ok=True)
                for old in case['saved_input_paths'].values():
                    path=Path(old);path.parent.mkdir(parents=True,exist_ok=True)
                    nib.save(image,path);nib.save(image,path.with_name(path.name.replace('cope','varcope',1)))
                    parent=path.parent.parent;(parent/'design.lev').write_text('2\n')
                    nib.save(image,parent/'mask.nii.gz');nib.save(image,parent/'bg_image.nii.gz')
                    number=int(re.search(r'cope(\d+)',path.name)[1]);nib.save(image,parent/'stats'/f'tdof_t{number}.nii.gz')
            corr=json.loads((r.TR_RECORD/'corrected-analysis-input.json').read_text());corr=json.loads(translate(json.dumps(corr)))
            candidate=root/'corrected/corrected_TR_1p710.feat';(candidate/'stats').mkdir(parents=True)
            original=Path(corr['original_feat'])
            candidate_fsf=r.replace_settings((original/'design.fsf').read_text(),{'fmri(tr)':'1.710'})
            (candidate/'design.fsf').write_text(candidate_fsf)
            for name in ('design.con','confoundevs.txt'):shutil.copyfile(original/name,candidate/name)
            for cope in (11,12):
                nib.save(image,candidate/'stats'/f'cope{cope}.nii.gz');nib.save(image,candidate/'stats'/f'varcope{cope}.nii.gz')
                nib.save(image,candidate/'stats'/f'tdof_t{cope}.nii.gz');nib.save(image,candidate/'mask.nii.gz');nib.save(image,candidate/'mean_func.nii.gz')
                for c in translated:
                    old=next(v for v in c['saved_input_paths'].values() if '/sub-11031/' in v)
                    cnum=int(re.search(r'-cope-(\d+)_',Path(c['path']).name)[1])
                    with patch.object(r,'ROOT',root):new,_=r.choose_input(old,cnum,candidate)
                    new.parent.mkdir(parents=True,exist_ok=True);nib.save(image,new);nib.save(image,new.with_name(new.name.replace('cope','varcope',1)))
                    nib.save(image,new.parent.parent/'mask.nii.gz');nib.save(image,new.parent.parent/'mean_func.nii.gz')
                    nib.save(image,new.parent/f'tdof_t{cnum}.nii.gz');(new.parent.parent/'design.lev').write_text('2\n')
            corr['active_feat']=str(candidate);(tr/'corrected-analysis-input.json').write_text(json.dumps(corr))
            with patch.object(r,'ROOT',root),patch.object(r,'TR_RECORD',tr),patch.dict(os.environ,{'FSLDIR':str(fsldir)}):
                planned=r.plan(ar,out,rec);again=r.plan(ar,out,rec)
            self.assertEqual(planned,again);self.assertEqual(len(planned),6)
            for case in planned:
                inputs=case['signature']['inputs'];self.assertEqual(len(inputs),132)
                self.assertEqual(sum(v['reason']!='reuse_original_participant_estimate' for v in inputs),2)
                old=r.core.fsl_settings(rec/'sources'/(case['name']+'-original.fsf'));new=r.core.fsl_settings(Path(case['fsf']))
                changed={k for k in old if old[k]!=new[k]}
                self.assertEqual(len([k for k in changed if k.startswith('feat_files(')]),132)
                self.assertFalse(any(k.startswith('fmri(evg') or k.startswith('fmri(con') for k in changed))
                self.assertTrue(changed<= {'fmri(outputdir)','fmri(regstandard)',*[k for k in old if k.startswith('feat_files(')]})

    def test_staging_expands_original_DOF_without_modifying_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'original.feat';(source/'stats').mkdir(parents=True);(source/'reg').mkdir()
            data=np.arange(60,dtype='float32').reshape(3,4,5);image=nib.Nifti1Image(data,np.eye(4))
            for name in ('stats/cope11.nii.gz','stats/varcope11.nii.gz','mask.nii.gz','mean_func.nii.gz'):nib.save(image,source/name)
            (source/'stats/dof').write_text('72\n');np.savetxt(source/'reg/example_func2standard.mat',np.eye(4))
            before={str(p):r.core.identity(p) for p in source.rglob('*') if p.is_file()}
            target=source/'stats/cope11.nii.gz';out=root/'staged.feat'
            a,receipt=r.stage_input(target,out,'reuse_original_participant_estimate',image)
            self.assertEqual(a,out/'stats/cope11.nii.gz');self.assertEqual((out/'design.lev').read_text(),'2\n')
            np.testing.assert_allclose(nib.load(out/'stats/tdof_t11.nii.gz').get_fdata(),72)
            self.assertEqual(r.core.sha(a),r.core.sha(target));self.assertEqual(before,{str(p):r.core.identity(p) for p in source.rglob('*') if p.is_file()})
            self.assertEqual(r.stage_input(target,out,'reuse_original_participant_estimate',image)[1],receipt)

    def test_group_completion_and_resume_validate_outputs_and_preserve_partial(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'old.gfeat';folder=source/'cope1.feat';(folder/'stats').mkdir(parents=True)
            for name,x in [('design.mat',np.ones((6,1))),('design.con',np.array([[1.],[-1.]])),('design.grp',np.ones((6,1)))]:
                with (folder/name).open('w') as stream:
                    stream.write('/Matrix\n');np.savetxt(stream,x)
            data=np.arange(60,dtype='float32').reshape(3,4,5)
            nib.save(nib.Nifti1Image(np.ones_like(data),np.eye(4)),folder/'mask.nii.gz')
            for sign in (1,2):
                for name in (f'stats/zstat{sign}.nii.gz',f'cluster_mask_zstat{sign}.nii.gz',f'thresh_zstat{sign}.nii.gz'):
                    nib.save(nib.Nifti1Image(data,np.eye(4)),folder/name)
            record=root/'record';record.mkdir();(record/'group_maps').mkdir()
            target=root/'new.gfeat';target.mkdir();(target/'partial.txt').write_text('retain')
            case=dict(name='ones-cope-11',original_group=str(source),output=str(target),fsf='/fixture.fsf',signature=dict(policy=r.POLICY,inputs=[]))
            def fake_command(argv,logdir):shutil.copytree(source,target)
            with patch.object(r.core,'run',side_effect=fake_command) as command:
                self.assertFalse(r.run_case(case,record,Path('/fake/feat'))['reused_completed_group'])
                self.assertTrue(r.run_case(case,record,Path('/fake/feat'))['reused_completed_group'])
                self.assertEqual(command.call_count,1)
                self.assertEqual(len(list(root.glob('new.gfeat.interrupted-*/partial.txt'))),1)
                (target/'cope1.feat/stats/zstat1.nii.gz').write_bytes(b'changed')
                with self.assertRaises(ValueError):r.run_case(case,record,Path('/fake/feat'))


if __name__=='__main__':unittest.main()
