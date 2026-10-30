import argparse
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

import audit as a


class SeedAuditTests(unittest.TestCase):
    def test_model_column_agreement_and_unsupported_filtering(self):
        rng=np.random.default_rng(3);x=14000+rng.normal(size=90);dm=rng.normal(size=(90,16));dm[:,11]=(x-x.mean())*.3
        settings=a.core.fsl_settings(a.ROOT/'templates/L1_task-ugr_model-3_type-ppi.fsf')
        self.assertEqual(a.check_model(settings,12,x,dm)['status'],'matches_up_to_positive_centering_scaling')
        self.assertEqual(a.check_model(settings,12,x[::-1],dm)['status'],'does_not_match_supported_model')
        self.assertEqual(a.check_model(dict(settings,**{'fmri(temphp_yn)':'1'}),12,x,dm)['status'],'not_checked_unsupported_processing')

    def test_raw_match_is_stricter_than_affine_agreement(self):
        x=14000+np.random.default_rng(4).normal(size=90)
        self.assertTrue(a.raw_match(x,x+1e-7)['raw_match'])
        result=a.raw_match(x,2*x+15)
        self.assertEqual(result['status'],'matches_up_to_positive_centering_scaling');self.assertFalse(result['raw_match'])

    def test_exact_field_of_view_embedding_has_no_interpolation(self):
        transform=np.eye(4);transform[:3,3]=[2,3,1]
        mask=dict(image=nib.Nifti1Image(np.zeros((3,4,5)),transform),indices=np.array([[0,0,0],[1,2,3]]))
        target=nib.Nifti1Image(np.zeros((10,10,10,20)),np.eye(4))
        np.testing.assert_array_equal(np.column_stack(a.exact_indices(mask,target)),[[2,3,1],[3,5,4]])
        transform[0,3]=2.5;mask['image']=nib.Nifti1Image(np.zeros((3,4,5)),transform)
        with self.assertRaises(ValueError):a.exact_indices(mask,target)

    def test_archive_preserves_case_distinct_sources_and_deduplicates_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);arc=a.Archive(root)
            first=arc.add_bytes(b'first','seed-AIns.txt','fixture');second=arc.add_bytes(b'second','seed-aIns.txt','fixture');third=arc.add_bytes(b'first','another.txt','fixture');arc.close()
            self.assertNotEqual(first['member'],second['member']);self.assertEqual(first['member'],third['member'])
            with tarfile.open(root/first['archive']) as tar:
                self.assertEqual(tar.extractfile(first['member']).read(),b'first');self.assertEqual(len(tar.getmembers()),2)

    def test_end_to_end_original_chain_and_stale_physiological_file(self):
        real=a.ROOT
        templates={task:a.core.fsl_settings(real/path) for task,path in [('ugr','templates/L1_task-ugr_model-3_type-ppi.fsf'),('trust','templates/archive/L1_task-trust_model-01_type-ppi.fsf')]}
        registry=json.loads(a.REGISTRY.read_text());rng=np.random.default_rng(9);shape=(6,7,8);t=60
        data=(14000+rng.normal(size=(*shape,t))*50).astype('float32');ijk={'AIns':(2,3,4),'AIns-vt':(3,4,5)}
        def write_settings(path,settings):
            path.write_text('\n'.join(f'set {k} "{v}"' for k,v in settings.items())+'\n')
        def matrix(path,values):
            with path.open('w') as stream:stream.write('/Matrix\n');np.savetxt(stream,values)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'masks').mkdir();(root/'code').mkdir();record=root/'record';record.mkdir()
            (root/'code/L1stats-ugr.sh').write_text('fslmeants -i /original/BOLD -m /original/seed-AIns.nii.gz -o /original/PHYS.txt\n')
            subprocess.run(['git','init','-q',str(root)],check=True);subprocess.run(['git','add','code/L1stats-ugr.sh'],cwd=root,check=True)
            subprocess.run(['git','-c','user.name=Fixture','-c','user.email=fixture@example.invalid','commit','-qm','fixture'],cwd=root,check=True)
            reg=root/'registry.json';reg.write_text(json.dumps(registry));bold=root/'original-BOLD.nii.gz';nib.save(nib.Nifti1Image(data,np.eye(4)),bold)
            session=root/'derivatives/fsl/sub-10369/ses-01';session.mkdir(parents=True)
            directories={}
            for seed,index in ijk.items():
                mask=np.zeros(shape,dtype='uint8');mask[index]=1;nib.save(nib.Nifti1Image(mask,np.eye(4)),root/'masks'/f'seed-{seed}.nii.gz')
            for task in ('ugr','trust'):
                for seed,index in ijk.items():
                    directory=session/f'L1_task-{task}_model-3_type-ppi_seed-{seed}_run-1_sm-5.feat';directory.mkdir();directories[task,seed]=directory
                    settings=dict(templates[task]);j=next(int(k[len('fmri(evtitle'):-1]) for k,v in settings.items() if k.startswith('fmri(evtitle') and v.lower()=='phys')
                    tc=session/f'sub-10369_ses-01_task-{task}_run-1_{seed}.txt';x=data[index].astype(float)
                    # One live external file is stale/overwritten; its saved model
                    # still contains the original B-series, so matching raw A must
                    # never be reported as historical B-seed confirmation.
                    current=data[ijk['AIns']].astype(float) if (task,seed)==('trust','AIns-vt') else x
                    np.savetxt(tc,current,fmt='%.12g');settings.update({'feat_files(1)':str(bold),f'fmri(custom{j})':str(tc),'fmri(npts)':str(t)})
                    write_settings(directory/'design.fsf',settings);dm=rng.normal(size=(t,60));dm[:,j-1]=(x-x.mean())*.3;matrix(directory/'design.mat',dm)
            second=session/'L2_task-ugr_model-3_type-ppi_seed-AIns_sm-5.gfeat';second.mkdir();write_settings(second/'design.fsf',{'feat_files(1)':str(directories['ugr','AIns'])})
            # A misleading activation-labelled group must still be followed
            # when its saved inputs point to the insula PPI analyses.
            group=root/'derivatives/fsl/L3-act/L3_type-act_n132_group-AIns/main.gfeat';group.mkdir(parents=True)
            write_settings(group/'design.fsf',{'feat_files(1)':str(second/'cope12.feat/stats/cope1.nii.gz')})
            state={}
            with patch.object(a,'ROOT',root),patch.object(a,'REGISTRY',reg):a.run(record,argparse.Namespace(max_bolds=2),state)
            self.assertEqual(state['physiological_EVs'],4)
            self.assertEqual(state['model_column_status']['does_not_match_supported_model'],1)
            self.assertEqual(state['numerical_reproduction']['model_supported_matches'],3)
            designs=json.loads((record/'saved_designs.json').read_text())
            target=next(row for row in designs if row.get('path')==str(directories['ugr','AIns']/'design.fsf'))
            self.assertTrue(target['referenced_by_discovered_group_chain'])
            comparisons=json.loads((record/'candidate_timecourse_comparisons.json').read_text())
            stale=[row for row in comparisons if row.get('task')=='trust' and row.get('folder_seed_label')=='AIns-vt' and row.get('raw_match')]
            self.assertEqual(len(stale),1);self.assertFalse(stale[0]['supports_saved_model_and_current_original_BOLD'])
            sources=json.loads((record/'source_manifest.json').read_text())
            self.assertFalse(any(row['source']==str(bold) for row in sources))
            self.assertTrue(any(row['role']=='historical_script' for row in sources))
            # Recover the actual original series from FEAT's snapshot even when
            # its external source was overwritten after the model was fitted.
            stale_dir=directories['trust','AIns-vt'];settings=a.core.fsl_settings(stale_dir/'design.fsf')
            j=next(int(k[len('fmri(evtitle'):-1]) for k,v in settings.items() if k.startswith('fmri(evtitle') and v.lower()=='phys')
            (stale_dir/'custom_timing_files').mkdir();np.savetxt(stale_dir/'custom_timing_files'/f'ev{j}.txt',data[ijk['AIns-vt']].astype(float),fmt='%.12g')
            rescued=root/'rescued';rescued.mkdir();state2={}
            with patch.object(a,'ROOT',root),patch.object(a,'REGISTRY',reg):a.run(rescued,argparse.Namespace(max_bolds=2),state2)
            self.assertEqual(state2['model_column_status']['matches_up_to_positive_centering_scaling'],4)
            self.assertEqual(state2['numerical_reproduction']['model_supported_matches'],4)
            self.assertEqual(state2['external_vs_cache_mismatches'],1)
            rows=json.loads((rescued/'physiological_EVs.json').read_text());cached=next(v for v in rows if v['task']=='trust' and v['folder_seed_label']=='AIns-vt')
            self.assertEqual(cached['series_source'],'FEAT_cached_EV');self.assertEqual(cached['external_model_check']['status'],'does_not_match_supported_model')


if __name__=='__main__':unittest.main()
