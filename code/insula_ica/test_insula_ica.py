"""Numerical safety tests and optional real-FSL residual-route integration."""
from pathlib import Path
import csv
import json
import os
import sys
import tempfile
import unittest
import numpy as np
import nibabel as nib
import insula_ica as p
import feat_residuals as f
import audit_linux2 as startup


def vest(path,matrix):
    matrix=np.atleast_2d(matrix)
    with path.open('w') as stream:
        stream.write(f'/NumWaves {matrix.shape[1]}\n/NumPoints {matrix.shape[0]}\n/Matrix\n')
        np.savetxt(stream,matrix,fmt='%.12g')


class NumericalTests(unittest.TestCase):
    def test_output_only_fsf(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'source.fsf';original='set fmri(outputdir) "/old/run.feat"\nset fmri(smooth) 5\nset fmri(prewhiten_yn) 1\n'
            path.write_text(original)
            revised=f.output_only_fsf(path,Path('/new/run.feat'))
            self.assertEqual(revised.replace('/new/run.feat','/old/run.feat'),original)
            self.assertEqual(path.read_text(),original)

    def test_multivariate_task_projection(self):
        rng=np.random.default_rng(2);maps=rng.normal(size=(130,5));maps[:,1]+=.8*maps[:,0]
        design=np.column_stack([np.ones(130),maps-maps.mean(axis=0)])
        weights=np.linalg.pinv(design)[1:];effect=np.array([1,-2,3,4,-5])
        cope=design[:,1:]@effect+15
        np.testing.assert_allclose(weights@cope,effect,atol=1e-12)
        # Spatial mixing commutes with a fixed task contrast, without any new GLM.
        beta=rng.normal(size=(130,7));contrast=rng.normal(size=7)
        np.testing.assert_allclose(weights@(beta@contrast),(weights@beta)@contrast,atol=1e-12)

    def test_partial_connectivity_conditions_all_components(self):
        rng=np.random.default_rng(3);tc=rng.normal(size=(150,5));tc[:,1]+=.7*tc[:,0]
        true=rng.normal(size=(5,8));target=tc@true
        np.testing.assert_allclose(np.linalg.pinv(tc)@target,true,atol=1e-12)

    def test_orientation(self):
        data=np.zeros((2,2,2,2));data[0,0,0]=[-5,4];data[1,1,1]=[2,-3]
        _,signs=p.orient_maps(data,np.ones((2,2,2),bool));self.assertEqual(signs,[-1,1])

    def test_qc_rejects_outdated_137_table(self):
        row={'sub_id':'10001','include_final':'TRUE',**{k:'TRUE' for k in startup.FLAGS}}
        with self.assertRaisesRegex(ValueError,'Corrected'):
            startup.cohort_view([row,dict(row,sub_id='10002')],{'10001'})
        self.assertEqual(set(startup.cohort_view([row,dict(row,sub_id='10002',include_final='FALSE')],{'10001'})),{'10001'})

    def test_no_silent_raw_route(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'config.json';path.write_text('{}')
            with self.assertRaisesRegex(ValueError,'retired'):p.config(path)

    def test_replay_changes_output_only(self):
        with tempfile.TemporaryDirectory() as folder:
            old=Path(folder)/'old.feat';old.mkdir();new=Path(folder)/'new.feat'
            text='film_gls --in=filtered_func_data --rn=stats --pd=design.mat --thr=1000 --sa --ms=5 --con=design.con'
            (old/'report_log.html').write_text('<pre>'+text+'</pre>')
            cmd,_=f.replay_command(old,new,True)
            self.assertIn('--rn='+str(new/'stats'),cmd)
            self.assertIn('--sa',cmd);self.assertFalse(new.exists())
            self.assertEqual((old/'report_log.html').read_text(),'<pre>'+text+'</pre>')

    def test_safe_recovery_routes(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);r,b,c,e,y=[root/name for name in ('residual','bold','conf','ev','filtered')]
            self.assertEqual(f.recovery_route(r,b,c,[e],y,['film_gls']),'blocked')
            y.touch();self.assertEqual(f.recovery_route(r,b,c,[e],y,['film_gls']),'replay_film')
            for path in (b,c,e):path.touch()
            self.assertEqual(f.recovery_route(r,b,c,[e],y,['film_gls']),'full_feat')
            r.touch();self.assertEqual(f.recovery_route(r,b,c,[e],y,None),'reuse')

    def test_historical_output_root_rejected(self):
        for bad in (p.ROOT/'derivatives/fsl',p.ROOT/'derivatives/fsl/sub-10001/anything',p.ROOT/'new.feat/child'):
            with self.assertRaisesRegex(ValueError,'cannot be inside'):
                f.safe_roots({'output_root':str(bad),'derivatives_root':str(p.ROOT/'derivatives/insula_ica/new')})

    def test_family(self):
        np.testing.assert_allclose(p.bonferroni([.001,.1]),[.03,1])


def film_probe(folder):
    folder.mkdir()
    rng=np.random.default_rng(11);t=180;tr=1.615
    dm=np.column_stack([np.sin(np.arange(t)/12),np.cos(np.arange(t)/22)])
    dm-=dm.mean(axis=0)
    eps=rng.normal(size=t)
    for i in range(1,t):eps[i]+=.8*eps[i-1]
    y=1000+dm@np.array([3.,-1.])+eps
    data=np.broadcast_to(y,(5,5,5,t)).copy().astype('float32')
    affine=np.diag([2.,2.,2.,1.]);im=nib.Nifti1Image(data,affine);im.header.set_zooms((2,2,2,tr));im.header.set_xyzt_units('mm','sec')
    nib.save(im,folder/'filtered_func_data.nii.gz');nib.save(nib.Nifti1Image(np.ones((5,5,5),'uint8'),affine),folder/'mask.nii.gz')
    vest(folder/'design.mat',dm);vest(folder/'design.con',np.eye(2))
    p.run(['film_gls','--in='+str(folder/'filtered_func_data.nii.gz'),'--rn='+str(folder/'stats'),'--pd='+str(folder/'design.mat'),
           '--con='+str(folder/'design.con'),'--thr=0','--outputPWdata'],folder/'commands')
    config={'output_root':str(folder),'residual_projection_tolerance':1e-4}
    row={'key':'synthetic_pw_probe','tr':tr,'prewhiten':'1','filtered_data':str(folder/'filtered_func_data.nii.gz')}
    result=f.residual_qc(folder,folder,row,config)
    assert result['scale']=='voxelwise_prewhitened'
    assert result['unwhitened_joint_fit_vs_res4d_relative_error']>.01
    assert result['max_whitened_decomposition_relative_error']<1e-4
    p.write_json(folder/'probe_result.json',result)
    print('FILM SCALE PROBE PASSED',result,flush=True)


def integration(folder):
    p.require(not folder.exists(),'Use a new synthetic output directory')
    original_root=p.ROOT;folder.mkdir();film_probe(folder/'film_scale_probe')
    p.ROOT=folder;startup.ROOT=folder
    import subprocess
    subprocess.run(['git','init','-q',str(folder)],check=True)
    subprocess.run(['git','-C',str(folder),'-c','user.name=Synthetic','-c','user.email=test@example.invalid','commit','--allow-empty','-qm','Synthetic fixture'],check=True)
    (folder/'code').mkdir();(folder/'code/sublist_n132.txt').write_text('10001\n10002\n10003\n')
    rng=np.random.default_rng(2026);shape=(49,47,37);t=180;tr=1.615
    affine=np.diag([2.,2.,2.,1.]);affine[:3,3]=[-48,-40,-24]
    ref=nib.Nifti1Image(np.zeros(shape,'uint8'),affine)
    atlas=Path(os.environ['FSLDIR'])/'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'
    xml=Path(os.environ['FSLDIR'])/'data/atlases/HarvardOxford-Cortical.xml'
    masks,_,_=p.anatomical_masks(atlas,xml,ref);brain=masks['bilateral'].copy()
    locations=rng.choice(brain.size,600,replace=False);brain.flat[locations]=True
    spatial=np.zeros((*shape,9),dtype='float32');labels=rng.choice(9,size=shape,p=np.arange(9,0,-1)/45)
    for k in range(9):spatial[...,k]=(labels==k)*rng.choice([-1.,1.],size=shape)
    records=[]
    for sub in ('10001','10002','10003'):
        qc={'sub_id':sub,'include_final':'TRUE'}
        for task,(model,_,_) in p.TASKS.items():
            qc[f'{task}_run1_full_ok']='TRUE';qc[f'{task}_run2_full_ok']='FALSE'
            key=f'sub-{sub}_ses-01_task-{task}_run-1'
            old=folder/f'derivatives/fsl/sub-{sub}/ses-01/L1_task-{task}_model-{model}_type-act_run-1_sm-5.feat';old.mkdir(parents=True)
            settings=(original_root/f'templates/L1_task-{task}_model-{model}_type-act.fsf').read_text()
            settings=settings.replace('NVOLUMES',str(t)).replace('OUTPUT',str(old)).replace('CONFOUNDEVS',str(old/'confoundevs.txt')).replace('DATA',str(folder/'missing_fmriprep24'/f'{key}.nii.gz')).replace('SMOOTH','5').replace('SHAPE_EV','10').replace('MISSED_TRIAL',str(old/'absent_miss.txt'))
            settings=__import__('re').sub(r'^set fmri\(prewhiten_yn\).*$', 'set fmri(prewhiten_yn) 0',settings,flags=__import__('re').M)
            nev=10 if task=='trust' else 11
            for j in range(1,nev):
                ev=old/f'ev{j}.txt';onsets=np.sort(rng.uniform(1,t*tr-20,12))
                np.savetxt(ev,np.column_stack([onsets,np.full(12,2),rng.uniform(.5,2,12)]))
                settings=__import__('re').sub(r'^set fmri\(custom'+str(j)+r'\).*$',f'set fmri(custom{j}) "{ev}"',settings,flags=__import__('re').M)
            (old/'design.fsf').write_text(settings)
            conf=rng.normal(size=(t,2));np.savetxt(old/'confoundevs.txt',conf)
            p.run(['feat_model',old/'design',old/'confoundevs.txt'],old/'fixture_commands')
            dm=p.fsl_matrix(old/'design.mat')
            latent=rng.laplace(size=(t,9)).astype('float32')
            signal=2*(spatial@latent.T)+rng.normal(size=(*shape,t))
            y=1000+signal+np.einsum('tp,vp->vt',dm,rng.normal(size=(brain.sum(),dm.shape[1]))).mean(axis=0)[None,None,None,:]
            data=np.zeros((*shape,t),dtype='float32');data[brain]=y[brain]
            image=nib.Nifti1Image(data,affine);image.header.set_zooms((2,2,2,tr));image.header.set_xyzt_units('mm','sec')
            nib.save(image,old/'filtered_func_data.nii.gz');nib.save(nib.Nifti1Image(brain.astype('uint8'),affine),old/'mask.nii.gz')
            nib.save(nib.Nifti1Image(data.mean(axis=3),affine),old/'mean_func.nii.gz')
            command=['film_gls','--in=filtered_func_data','--rn=stats','--pd=design.mat','--con=design.con','--thr=0','--noest']
            # Run only in the synthetic historical folder to create the reference.
            result=subprocess.run(command,cwd=old,capture_output=True,text=True)
            p.require(result.returncode==0,'Fixture FILM failed: '+result.stderr)
            (old/'report_log.html').write_text('<pre>'+__import__('shlex').join(command)+'</pre>')
            (old/'stats/res4d.nii.gz').unlink()  # Simulate the historical cleanup.
        records.append(qc)
    qc_path=folder/'code/master_subject_table_full.csv'
    with qc_path.open('w') as stream:
        w=csv.DictWriter(stream,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    c=json.loads((original_root/'code/insula_ica/config.example.json').read_text())
    c.update(master_qc=str(qc_path),output_root=str(folder/'output'),derivatives_root=str(folder/'derivatives/insula_ica'),expected_cohort_size=3,
             ugr_pmods_increase_with_fairness=True,ugr_coding_evidence='SYNTHETIC amplitudes only')
    template=folder/'network.nii.gz';nib.save(nib.Nifti1Image(rng.normal(size=shape),affine),template)
    c['network_templates']=[{'name':'synthetic only','path':str(template),'citation':'synthetic fixture'}]
    before={str(x):p.sha(x) for x in (folder/'derivatives/fsl').rglob('*') if x.is_file()}
    assert p.audit(c)==0
    a=json.loads((folder/'output/audit/audit.json').read_text());assert a['recovery_counts']=={'replay_film':6}
    p.reconstruct(c,'pilot',False)
    assert not (Path(c['derivatives_root'])/'feat_reconstruction').exists()
    p.reconstruct(c,'pilot',True);assert p.verify(c,'pilot')==0
    p.reconstruct(c,'full',True);assert p.verify(c,'full')==0
    for path,digest in before.items():assert p.sha(path)==digest,'Historical result modified'
    p.prepare(c);p.ica(c,'pilot');p.regress(c,'pilot');p.summarize(c,'pilot')
    assert len(p.read_table(folder/'output/pilot/group_task_effects.tsv'))==30
    p.ica(c,'full');p.regress(c,'full');p.summarize(c,'full')
    # A reviewed prewhitened fit is deliberately required, never inferred from the file name.
    c['residual_scale_review']['accept_voxelwise_prewhitened']=False
    print('FEAT RESIDUAL INTEGRATION PASSED:',folder/'output/full/results_summary.md',flush=True)


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--integration':integration(Path(sys.argv[2]).resolve())
    else:unittest.main()
