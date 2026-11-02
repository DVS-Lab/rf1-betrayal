#!/usr/bin/env python3
"""Masked stage-1 / whole-brain stage-2 dual regression on completed 3/5 ICAs."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import fcntl
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import traceback
import uuid
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
from scipy import ndimage

import run as f
import network_summary as summary

POLICY='masked-functional-DR-wholebrain-3-5-sm0-sm5-paired-scales-v1'
SUCCESS='wholebrain_dual_regression_complete_review_pending'
FAMILY_TESTS=64  # 32 prespecified components, positive and negative tails.


def vest(path,values,contrast=False):
    values=np.atleast_2d(values)
    if path.is_file():
        old=f.g.core.fsl_matrix(path);f.g.core.require(old.shape==values.shape and np.array_equal(old,values),'Frozen group design/contrast changed')
        return
    with path.open('w') as stream:
        stream.write(f'/NumWaves {values.shape[1]}\n'+('/NumContrasts ' if contrast else '/NumPoints ')+str(values.shape[0])+'\n/Matrix\n')
        np.savetxt(stream,values,fmt='%.12g')


def temporal_operator(weights,residual):
    tc=(weights@np.asarray(residual,dtype=float)).T;tc-=tc.mean(axis=0)
    sd=tc.std(axis=0,ddof=1)
    f.g.core.require(np.isfinite(tc).all() and (sd>1e-8).all(),'Invalid masked residual timecourse')
    z=tc/sd;condition=float(np.linalg.cond(z))
    f.g.core.require(np.linalg.matrix_rank(z)==z.shape[1] and condition<1e6,'Nonestimable component timecourses')
    return np.linalg.pinv(z),z,dict(condition=condition,SD_before_normalization=sd.tolist(),TC_mean_after_normalization=z.mean(axis=0).tolist())


def bold_projection(operator,filtered,design,pes):
    """Use existing joint FILM PEs, without refitting task/confounds."""
    y=np.asarray(filtered,dtype=float);pe=np.asarray(pes,dtype=float)
    return operator@(y-y.mean(axis=1,keepdims=True)).T-(operator@design)@pe.T


def validate_run(arm,row,source_root):
    marker=source_root/'run_arrays'/(row['key']+'_complete.json');cp=f.read(marker)
    f.g.core.require(cp['signature']['mask_sha256']==f.MASK_SHA,'Residual cache mask changed')
    f.g.core.require(all(f.g.core.identity(p['path'])==p for p in cp['arrays']),'Residual cache changed')
    f.g.core.require(all(f.g.core.identity(p['path'],p['sha256'] is not None)==p for p in cp['qc']['source_identities']),'Retained fitted products changed')
    return dict(arm=arm,key=row['key'],residual_checkpoint_sha256=f.g.core.sha(marker),sources=cp['qc']['source_identities'])


def analyze_run(row,source_root,folder,families,original,ica_mask,wb,reference):
    """One whole-brain read per representation for BOTH requested dimensions."""
    feat=Path(row['active_feat']);design=f.g.core.fsl_matrix(feat/'design.mat');diagnostics=[]
    for mode in f.g.inputs.MODES:
        selected={k:v for k,v in families.items() if v['mode']==mode};operators={}
        residual=np.load(source_root/'run_arrays'/(row['key']+'_'+mode+'.npy'))[ica_mask[original]]
        for key,item in selected.items():
            op,z,qc=temporal_operator(item['weights'],residual);operators[key]=op
            np.savetxt(folder/(key+'_timecourses.tsv'),z,delimiter='\t',fmt='%.12g');diagnostics.append(dict(family=key,**qc))
        matrices={k:np.empty((int(wb.sum()),v['dim']),dtype='float32') for k,v in selected.items()}
        if mode=='native_film':
            image=nib.load(row['native_residual']);data=np.asarray(image.dataobj,dtype='float32')[wb]
            f.g.core.require(data.shape[1]==row['nvol'] and np.isfinite(data).all(),'Nonfinite/mistimed native whole-brain residuals')
            data=data-data.mean(axis=1,keepdims=True)
            for key,operator in operators.items():matrices[key][:]=(operator@data.T).T
            del data
        else:
            data=np.asarray(nib.load(row['filtered_data']).dataobj,dtype='float32')[wb]
            pes=np.column_stack([np.asarray(nib.load(feat/f'stats/pe{k+1}.nii.gz').dataobj,dtype='float32')[wb] for k in range(design.shape[1])])
            f.g.core.require(np.isfinite(data).all() and np.isfinite(pes).all(),'Nonfinite whole-brain fitted values/PEs')
            for start in range(0,len(data),4096):
                for key,operator in operators.items():matrices[key][start:start+4096]=bold_projection(operator,data[start:start+4096],design,pes[start:start+4096]).T
            del data,pes
        for key,values in matrices.items():
            f.g.core.require(np.isfinite(values).all(),'Nonfinite dual-regression coefficient')
            full=np.zeros((*wb.shape,values.shape[1]),dtype='float32');full[wb]=values
            f.g.core.save_image(folder/(key+'_DR.nii.gz'),full,reference)
    return dict(diagnostics=diagnostics,all_components_fitted_together=True,stage1_from_masked_matching_residuals=True,new_nuisance_or_task_fit=False,new_temporal_filter=False)


def randomise_command(merged,prefix,mask,design,contrasts,permutations,seed=20261101):
    return ['randomise','-i',str(merged),'-o',str(prefix),'-m',str(mask),'-d',str(design),'-t',str(contrasts),'-n',str(permutations),'-T','--seed='+str(seed)]


def randomise_job(merged,folder,mask,design,contrasts,permutations,components):
    command=randomise_command(merged,folder/'dr',mask,design,contrasts,permutations)
    signature=dict(command=command,merged=f.g.core.identity(merged),mask=f.g.core.identity(mask),design=f.g.core.identity(design),contrasts=f.g.core.identity(contrasts),family_tests=FAMILY_TESTS,components=components)
    marker=folder/'complete.json'
    if marker.is_file():
        old=f.read(marker);f.g.core.require(old['signature']==signature and all(f.g.core.identity(v['path'])==v for v in old['outputs']),'Completed randomise inputs/outputs changed')
        return dict(state='complete',reused=True,folder=str(folder))
    if folder.exists():
        commands=subprocess.check_output(['ps','-u',str(os.getuid()),'-o','args='],text=True)
        f.g.core.require(not any(str(folder/'dr') in line for line in commands.splitlines()),'An earlier randomise process is still active for this output')
        f.preserve(folder)
    folder.mkdir(parents=True);f.save(folder/'command.json',signature)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (folder/'randomise.log').open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    outputs=[folder/f'dr_{kind}_tstat{tail}.nii.gz' for tail in (1,2) for kind in ('tfce_corrp',)]
    outputs.extend(folder/f'dr_tstat{tail}.nii.gz' for tail in (1,2))
    f.g.core.require(all(p.is_file() for p in outputs),'Missing positive/negative randomise outputs')
    for p in outputs:
        im=nib.load(p);data=im.get_fdata();f.g.core.require(im.shape==nib.load(mask).shape and f.g.core.same_grid(im,nib.load(mask)) and np.isfinite(data).all(),'Invalid randomise geometry/values')
        if 'corrp' in p.name:f.g.core.require(data.min()>=0 and data.max()<=1,'Invalid 1-p output')
    f.save(marker,dict(signature=signature,outputs=[f.g.core.identity(p) for p in outputs]))
    return dict(state='complete',reused=False,folder=str(folder))


def common_mask(rows,reference,record,root,original,ica_mask):
    mask=np.ones(reference.shape,bool);union=np.zeros_like(mask)
    for arm,rr in rows.items():
        for row in rr:
            feat=Path(row['active_feat']);m=nib.load(feat/'mask.nii.gz');sigma=nib.load(feat/'stats/sigmasquareds.nii.gz')
            f.g.core.require(f.g.core.same_grid(m,reference) and f.g.core.same_grid(sigma,reference),'Whole-brain fitted geometry differs')
            values=sigma.get_fdata();valid=(m.get_fdata()>0)&np.isfinite(values)&(values>0);mask&=valid;union|=valid
    f.g.core.require(mask.sum()>1000 and not (ica_mask&~mask).any(),'Insufficient common whole-brain or defining-ROI support')
    infer=mask&~original;distance=ndimage.distance_transform_edt(~original,sampling=reference.header.get_zooms()[:3]);distant=infer&(distance>=10)
    for name,data in [('wholebrain_common_mask',mask),('inference_outside_original_ROI',infer),('distant_targets_10mm_mask',distant)]:
        f.g.freeze_nifti(root/(name+'.nii.gz'),data.astype('uint8'),reference.affine);f.d.copy_frozen(root/(name+'.nii.gz'),record/'group_maps'/(name+'.nii.gz'))
    f.save(record/'wholebrain_coverage.json',dict(common_voxels=int(mask.sum()),union_voxels=int(union.sum()),fraction_of_union=float(mask.sum()/union.sum()),inference_voxels=int(infer.sum()),distant_target_voxels=int(distant.sum()),definition='Common positive finite fitted sigma and FEAT mask over all 980 run-arm entries; all 4D values must subsequently be finite',defining_ROI_excluded_from_group_inference_and_network_correspondence=True,no_participant_or_voxel_selection_based_on_connectivity=True))
    return mask,infer,distant


def aggregate_family(item,rows,root,record,reference,wb,source_run_markers):
    key=item['key'];arm=item['arm'];rr=rows[arm];subjects=sorted({r['subject'] for r in rr});folder=root/'participants'/key;folder.mkdir(parents=True,exist_ok=True)
    for subject in subjects:
        marker=folder/(subject+'.complete.json');relevant=[r for r in rr if r['subject']==subject]
        signature=dict(run_receipts=[source_run_markers[arm,r['key']] for r in relevant],family=key,aggregation='equal runs within task; equal task means within participant')
        if marker.is_file():
            old=f.read(marker);f.g.core.require(old['signature']==signature and all(f.g.core.identity(p['path'])==p for p in old['outputs']),'Participant aggregation cache changed');continue
        taskmeans={}
        for task in ('ugr','trust'):
            entries=[r for r in relevant if r['task']==task];f.g.core.require(entries,'Missing participant/task')
            arrays=[nib.load(root/'runs'/arm/r['key']/(key+'_DR.nii.gz')).get_fdata(dtype='float32') for r in entries]
            taskmeans[task]=np.mean(arrays,axis=0,dtype=np.float64).astype('float32')
        outputs=[]
        for label,data in [*taskmeans.items(),('pooled',(taskmeans['ugr'].astype(float)+taskmeans['trust'])/2)]:
            path=folder/(subject+'_'+label+'.nii.gz');f.g.core.save_image(path,np.asarray(data,dtype='float32'),reference);outputs.append(f.g.core.identity(path))
        f.save(marker,dict(signature=signature,outputs=outputs))
    for label in ('pooled','ugr','trust'):
        target=record/'group_maps'/(key+'_'+label+'_mean.nii.gz')
        mean=np.zeros((*wb.shape,item['dim']),dtype='float64')
        for sub in subjects:mean+=nib.load(folder/(sub+'_'+label+'.nii.gz')).get_fdata(dtype='float32')/len(subjects)
        f.g.freeze_nifti(target,mean.astype('float32'),reference.affine)
    merges=root/'merged'/key;merges.mkdir(parents=True,exist_ok=True)
    for k in range(item['dim']):
        path=merges/f'IC{k+1}.nii.gz';marker=path.with_name(path.name+'.complete.json')
        signature=dict(subjects=subjects,participant_markers=[f.g.core.sha(folder/(s+'.complete.json')) for s in subjects],component=k+1)
        if marker.is_file():
            old=f.read(marker);f.g.core.require(old['signature']==signature and f.g.core.identity(path)==old['output'],'Merged participant input changed');continue
        values=np.zeros((*wb.shape,len(subjects)),dtype='float32')
        for index,subject in enumerate(subjects):values[...,index]=nib.load(folder/(subject+'_pooled.nii.gz')).get_fdata(dtype='float32')[...,k]
        f.g.core.save_image(path,values,reference);f.save(marker,dict(signature=signature,output=f.g.core.identity(path)))
    return subjects


def region_targets(item,mean,distant,reference,fsldir,networks,labels):
    rows=[]
    for name,xmlname in [('Cortical','cort'),('Subcortical','sub')]:
        path=fsldir/f'data/atlases/HarvardOxford/HarvardOxford-{xmlname}-maxprob-thr25-1mm.nii.gz'
        source=fsldir/f'data/atlases/HarvardOxford-{name}.xml'
        atlas=resample_from_to(nib.load(path),reference,order=0).get_fdata().astype(int)
        atlas_labels={int(x.attrib['index'])+1:x.text for x in ET.parse(source).findall('./data/label')}
        for index,label in atlas_labels.items():
            mask=(atlas==index)&distant
            if mask.sum()<20:continue
            if any(x in label.lower() for x in ('white matter','cerebral cortex','ventricle','brain-stem')):continue
            for k,value in enumerate(mean[mask].mean(axis=0)):rows.append(dict(family=item['key'],component=k+1,atlas=name,target=label,distant_voxels=int(mask.sum()),mean_background_coefficient=float(value),descriptive_only=True))
    corr=[]
    for k in range(item['dim']):
        y=mean[distant,k]
        for j,label in enumerate(labels):
            x=networks[distant,j]
            r=float(np.corrcoef(x,y)[0,1]) if np.std(x)>1e-10 and np.std(y)>1e-10 else None
            corr.append(dict(family=item['key'],component=k+1,Smith2009_RSN=j+1,reference_label=label,spatial_r=r,mask='Common brain, outside original ROI plus 10-mm guard band',interpretation='Descriptive spatial correspondence, not exclusive membership or inferential network assignment'))
    return rows,corr


def network_figure(item,mean,tmap,visible,wb,original,background,reference,path,label):
    dim=item['dim'];fig,axes=plt.subplots(dim,3,figsize=(12,2.65*dim),squeeze=False);vmax=np.percentile(background[background>0],99.5)
    pos=[reference.affine[k,3]+np.arange(reference.shape[k])*reference.affine[k,k] for k in range(3)]
    for k in range(dim):
        outside=wb&~original
        values=mean[...,k] if tmap is None else tmap[...,k]
        valid=outside if visible is None else outside&visible[...,k]
        candidates=np.where(valid,np.abs(values),-1);peak=np.unravel_index(np.argmax(candidates),values.shape)
        no_voxels=not valid.any()
        if no_voxels:peak=np.unravel_index(np.argmax(np.where(outside,np.abs(mean[...,k]),-1)),values.shape)
        cutoff=float(np.percentile(np.abs(values[outside]),70)) if visible is None else 0
        scale=max(float(np.percentile(np.abs(values[outside]),98)),1e-10)
        for axis,ax in enumerate(axes[k]):
            h,v=((1,2),(0,2),(0,1))[axis];xx,yy=pos[h],pos[v];extent=(xx[0],xx[-1],yy[0],yy[-1]);plane=lambda a:np.take(a,int(peak[axis]),axis=axis).T
            ax.imshow(plane(background),cmap='gray',origin='lower',extent=extent,vmin=0,vmax=vmax)
            selected=plane(valid)&(np.abs(plane(values))>cutoff)
            ax.imshow(np.ma.masked_where(~selected,plane(values)),cmap='coolwarm',vmin=-scale,vmax=scale,origin='lower',extent=extent,interpolation='nearest',alpha=.9)
            ax.set_title(f'IC{k+1}: {"xyz"[axis]}={pos[axis][peak[axis]]:.1f} mm'+(' (no surviving voxels)' if no_voxels else ''),fontsize=8);ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle(item['key']+' | '+label,fontsize=12)
    fig.text(.5,.01,'N=132, equal runs within task and equal task weights. Red positive, blue negative partial regression coefficient/statistic.\nDefining ROI excluded from display. '+('Descriptive: largest 30% absolute coefficients; not a significance threshold.' if visible is None else 'TFCE-corrected inference; signed t values shown only where the stated family threshold is passed.'),ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.95));fig.savefig(path,dpi=135);plt.close(fig)


def execute(args,record,root,state):
    sources,states,rows,reference,ica_mask,original,families=summary.load_sources(Path(state['sm0_record']),Path(state['sm5_record']))
    source_roots={arm:Path(s['output_root']) for arm,s in states.items()};fsldir=Path(os.environ['FSLDIR'])
    meta_path=Path(__file__).parent/'reference/Smith2009_rsn10_metadata.json';meta=f.read(meta_path)
    template=Path(__file__).parent/'reference/Smith2009_rsn10_analysis_grid.nii.gz'
    f.g.core.require(f.g.core.sha(template)==meta['derived_sha256'] and f.g.core.same_grid(nib.load(template),reference),'Independent network reference changed/grid differs')
    shared=dict(policy=POLICY,source_records={k:str(v) for k,v in sources.items()},registries={k:f.g.core.identity(p/'accepted-inputs.json') for k,p in sources.items()},maps={k:v['map_sha256'] for k,v in families.items()},permutations=args.permutations,family_tests=FAMILY_TESTS,
        independent_reference_sha256=f.g.core.sha(template),reference_metadata_sha256=f.g.core.sha(meta_path),
        implementation_sha256={Path(m.__file__).name:f.g.core.sha(m.__file__) for m in (sys.modules[__name__],summary,f,f.g,f.g.inputs)})
    f.d.freeze_json(record/'provenance.json',shared);software=f.g.core.software();software['executables']['randomise']=shutil.which('randomise')
    f.g.core.require(software['executables']['randomise'] is not None and software['fsl']==f.read(sources['sm0']/'software.json')['fsl'],'Randomise unavailable or FSL version differs from source ICA fits');f.save(record/'software.json',software)
    for key,item in families.items():f.d.copy_frozen(item['map_path'],record/'primary_maps'/(key+'.nii.gz'))
    f.d.copy_frozen(sources['sm0']/'primary_maps/ICA_common_functional_cluster.nii.gz',record/'primary_maps/ICA_mask.nii.gz')
    summary.anatomy(families,reference,ica_mask,record,fsldir);summary.task_summary(sources,record)
    state.update(state='anatomy_and_social_contrasts_ready',anatomy_ready=True);f.publish_stage(record,state,args.publish,'Publish functional component anatomy and paired social-specificity summaries')
    wb,infer,distant=common_mask(rows,reference,record,root,original,ica_mask)
    run_signature=dict(shared,wholebrain_mask_sha256=f.g.core.sha(root/'wholebrain_common_mask.nii.gz'))
    gate=f.g.inputs.MemoryGate();lock=threading.Lock();source_run_markers={}
    def worker(arm,row):
        key=arm+'/'+row['key'];folder=root/'runs'/arm/row['key'];marker=folder/'complete.json'
        try:
            source=validate_run(arm,row,source_roots[arm]);signature=dict(run_signature,run=source)
            if marker.is_file():
                old=f.read(marker);f.g.core.require(old['signature']==signature and all(f.g.core.identity(p['path'])==p for p in old['outputs']),'Completed whole-brain DR source/output changed');reused=True
            else:
                if folder.exists():f.preserve(folder)
                folder.mkdir(parents=True)
                with gate.reserve(4):qc=analyze_run(row,source_roots[arm],folder,{k:v for k,v in families.items() if v['arm']==arm},original,ica_mask,wb,reference)
                f.g.core.require(validate_run(arm,row,source_roots[arm])==source,'Source changed during whole-brain regression')
                outputs=[f.g.core.identity(p) for p in folder.iterdir() if p.is_file()];f.save(marker,dict(signature=signature,qc=qc,outputs=outputs));reused=False
            result=dict(state='complete',reused=reused);source_run_markers[arm,row['key']]=f.g.core.sha(marker)
            error=record/'errors'/(key.replace('/','_')+'.txt')
            if error.is_file():error.rename(record/'attempts'/(error.stem+'-'+uuid.uuid4().hex[:6]+'.txt'))
        except Exception as exc:result=dict(state='failed',error=str(exc));(record/'errors'/(key.replace('/','_')+'.txt')).write_text(traceback.format_exc())
        with lock:
            state['runs'][key]=result;f.save(record/'state.json',state)
            print(f'Dual regression: {sum(v["state"]=="complete" for v in state["runs"].values())}/980 | {key}: {result["state"]}',flush=True)
        return result
    state['state']='masked_to_wholebrain_dual_regression';f.save(record/'state.json',state)
    pilots=[rows['sm0'][0],next(r for r in rows['sm0'] if r['key']=='sub-11031_ses-01_task-ugr_run-1'),next(r for r in rows['sm0'] if r['key']=='sub-10817_ses-01_task-ugr_run-1')];pilot_keys={r['key'] for r in pilots}
    for arm,rr in rows.items():
        for row in rr:
            if row['key'] in pilot_keys:f.g.core.require(worker(arm,row)['state']=='complete','Whole-brain DR pilot failed; no expansion')
    state['pilot_passed']=True;f.save(record/'state.json',state)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(worker,arm,row) for arm,rr in rows.items() for row in rr if row['key'] not in pilot_keys]):future.result()
    f.g.core.require(len(state['runs'])==980 and all(v['state']=='complete' for v in state['runs'].values()),'Whole-brain runs incomplete; no incomplete-cohort group inference')
    state['state']='participant_aggregation_and_group_means';f.save(record/'state.json',state)
    subjects=None
    for key,item in families.items():subjects=aggregate_family(item,rows,root,record,reference,wb,source_run_markers);state['groups'][key]=dict(state='complete',participants=len(subjects));f.save(record/'state.json',state)
    vest(root/'group_design.mat',np.ones((len(subjects),1)));vest(root/'group_design.con',np.array([[1],[-1]]),contrast=True);f.save(record/'participant_order.json',subjects)
    f.save(record/'inference_policy.json',dict(permutations=args.permutations,TFCE=True,unit='132 participant maps, not 490 runs',maps=32,tails=2,all_tests=FAMILY_TESTS,per_component_tail_1_minus_p=.95,within_basis_1_minus_p={'dim3':1-.05/6,'dim5':1-.05/10},whole_requested_family_1_minus_p=1-.05/FAMILY_TESTS,primary_group_question='Pooled mean background partial coupling, equal task weights',per_task_means='Descriptive only, not social-specific connectivity',seed=20261101,common_mask=True,whole_original_127_voxel_ROI_excluded=True))
    f.publish_stage(record,state,args.publish,'Publish whole-brain dual-regression group means; start TFCE inference')
    state['state']='randomise_TFCE_group_networks';f.save(record/'state.json',state)
    def infer_component(key,item,k):
        job=key+f'/IC{k+1}'
        try:
            with gate.reserve(2):value=randomise_job(root/'merged'/key/f'IC{k+1}.nii.gz',root/'randomise'/key/f'IC{k+1}',root/'inference_outside_original_ROI.nii.gz',root/'group_design.mat',root/'group_design.con',args.permutations,item['dim'])
            error=record/'errors'/('randomise_'+job.replace('/','_')+'.txt')
            if error.is_file():error.rename(record/'attempts'/(error.stem+'-'+uuid.uuid4().hex[:6]+'.txt'))
        except Exception as exc:value=dict(state='failed',error=str(exc));(record/'errors'/('randomise_'+job.replace('/','_')+'.txt')).write_text(traceback.format_exc())
        with lock:state['randomise'][job]=value;f.save(record/'state.json',state);print('Randomise',job,value['state'],flush=True)
    with ThreadPoolExecutor(max_workers=args.randomise_jobs) as pool:
        for future in as_completed([pool.submit(infer_component,key,item,k) for key,item in families.items() for k in range(item['dim'])]):future.result()
    f.g.core.require(len(state['randomise'])==32 and all(v['state']=='complete' for v in state['randomise'].values()),'Randomise incomplete; completed maps retained')
    state['state']='summarizing_group_networks';f.save(record/'state.json',state)
    f.g.core.require(f.g.core.sha(template)==meta['derived_sha256'] and f.g.core.same_grid(nib.load(template),reference),'Independent network reference changed/grid differs')
    networks=nib.load(template).get_fdata();f.save(record/'independent_network_reference.json',meta)
    background=resample_from_to(nib.load(fsldir/'data/standard/MNI152_T1_2mm_brain.nii.gz'),reference,order=1).get_fdata()
    targets=[];correlations=[];counts=[]
    for key,item in families.items():
        mean=nib.load(record/'group_maps'/(key+'_pooled_mean.nii.gz')).get_fdata();stats_signed=np.zeros_like(mean);global_visible=np.zeros_like(mean,dtype=bool);within_visible=np.zeros_like(global_visible)
        for k in range(item['dim']):
            folder=root/'randomise'/key/f'IC{k+1}'
            for tail in (1,2):
                pmap=nib.load(folder/f'dr_tfce_corrp_tstat{tail}.nii.gz').get_fdata();tmap=nib.load(folder/f'dr_tstat{tail}.nii.gz').get_fdata();direction='positive' if tail==1 else 'negative'
                valid=infer&(tmap>0);primary=valid&(pmap>=1-.05/FAMILY_TESTS);basis=valid&(pmap>=1-.05/(2*item['dim']))
                stats_signed[...,k]+=np.where(valid,tmap*(1 if tail==1 else -1),0);global_visible[...,k]|=primary;within_visible[...,k]|=basis
                counts.append(dict(family=key,component=k+1,direction=direction,voxelwise_TFCE_per_map_FWE05_voxels=int((valid&(pmap>=.95)).sum()),within_basis_components_and_tails_FWE05_voxels=int(basis.sum()),all_32_components_and_two_tails_FWE05_voxels=int(primary.sum())))
                for name in (f'dr_tfce_corrp_tstat{tail}.nii.gz',f'dr_tstat{tail}.nii.gz'):f.d.copy_frozen(folder/name,record/'inference_maps'/f'{key}_IC{k+1}_{name}')
            f.d.copy_frozen(folder/'command.json',record/'commands'/f'{key}_IC{k+1}.json')
        for suffix,visible in [('all_models_FWE05',global_visible),('within_basis_FWE05',within_visible)]:
            f.g.freeze_nifti(record/'group_maps'/(key+'_'+suffix+'_signed_t.nii.gz'),np.where(visible,stats_signed,0).astype('float32'),reference.affine)
            network_figure(item,mean,stats_signed,visible,wb,original,background,reference,record/'figures'/(key+'_'+suffix+'_networks.png'),suffix)
        network_figure(item,mean,None,None,wb,original,background,reference,record/'figures'/(key+'_descriptive_mean_networks.png'),'descriptive mean coefficients')
        regional,corr=region_targets(item,mean,distant,reference,fsldir,networks,meta['labels']);targets+=regional;correlations+=corr
    summary.write_table(record/'inference_summary.tsv',counts);summary.write_table(record/'anatomical_network_targets.tsv',targets);summary.write_table(record/'independent_RSN_correspondence.tsv',correlations)
    state.update(state=SUCCESS,completed_components=32,completed_randomise_jobs=32);f.save(record/'state.json',state)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sm0-record',type=Path);p.add_argument('--sm5-record',type=Path);p.add_argument('--resume-record',type=Path);p.add_argument('--new-record',action='store_true');p.add_argument('--jobs',type=int,default=32);p.add_argument('--randomise-jobs',type=int,default=32);p.add_argument('--permutations',type=int,default=5000);p.add_argument('--publish',action='store_true');args=p.parse_args()
    f.g.core.require(1<=args.jobs<=64 and 1<=args.randomise_jobs<=64 and args.permutations>=5000,'Use 1–64 RAM-gated workers and at least 5000 permutations')
    args.sm0_record=args.sm0_record or sorted(f.RECORDS.glob('linux2-functional-cluster-*'))[-1];args.sm5_record=args.sm5_record or sorted(f.RECORDS.glob('linux2-functional-sm5-*'))[-1]
    if not args.resume_record and not args.new_record:
        for candidate in sorted(f.RECORDS.glob('linux2-functional-wholebrain-*'),reverse=True):
            if (candidate/'state.json').is_file():
                old=f.read(candidate/'state.json')
                if old.get('policy')==POLICY and Path(old['sm0_record']).resolve()==args.sm0_record.resolve() and Path(old['sm5_record']).resolve()==args.sm5_record.resolve() and old['permutations']==args.permutations:args.resume_record=candidate;break
    output=f.ROOT/'code/insula_ica/output';output.mkdir(exist_ok=True)
    with (output/'functional-wholebrain.lock').open('a') as lock,(output/'functional-cluster.lock').open('a') as a,(output/'functional-sm5.lock').open('a') as b:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:print('Whole-brain dual-regression runner already active.');return
        fcntl.flock(a,fcntl.LOCK_SH);fcntl.flock(b,fcntl.LOCK_SH)
        if args.resume_record:
            record=args.resume_record;state=f.read(record/'state.json');root=Path(state['output_root'])
            f.g.core.require(record.resolve().is_relative_to(f.RECORDS.resolve()) and root.resolve().is_relative_to((f.ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('functional-wholebrain-'),'Invalid network resume root')
        else:
            stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6];record=f.RECORDS/('linux2-functional-wholebrain-'+stamp);record.mkdir();root=f.ROOT/'derivatives/insula_ica'/('functional-wholebrain-'+stamp);root.mkdir(parents=True)
            state=dict(policy=POLICY,state='preflight',sm0_record=str(args.sm0_record),sm5_record=str(args.sm5_record),output_root=str(root),permutations=args.permutations,runs={},groups={},randomise={},pilot_passed=False)
        for name in ('figures','primary_maps','group_maps','inference_maps','commands','attempts','errors'):(record/name).mkdir(exist_ok=True)
        for name in ('runs','participants','merged','randomise'):(root/name).mkdir(exist_ok=True)
        if (record/'failure.txt').is_file():(record/'failure.txt').rename(record/'attempts'/('failure-'+uuid.uuid4().hex[:6]+'.txt'))
        state.pop('error',None);state.update(jobs=args.jobs,randomise_jobs=args.randomise_jobs,runner_pid=os.getpid(),code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=f.ROOT,text=True).strip());f.save(record/'state.json',state)
        (record/'report.md').write_text('# Whole-brain dual regression from functional masked ICA\n\nEight fixed bases: sm0/sm5 × BOLD-scale/native FILM × 3/5 components, 32 components in total. Stage 1 uses only the same 117 masked residual voxels, fits all maps plus a spatial intercept, then centers and variance-normalizes the component timecourses. Stage 2 regresses all timecourses for that basis together onto matching whole-brain residuals. No new nuisance/task fit, filter, whitening, ICA, FEAT, fMRIPrep or interpolation of participant data. Native FILM is an innovation-scale sensitivity, not ordinary unwhitened temporal covariance.\n\nRuns average within participant/task; tasks have equal weight in pooled participant maps. Randomise uses 132 independent maps, 5000 or more permutations, 3D TFCE and two signed contrasts. Defining 127-voxel ROI is excluded from inference. Report per-map TFCE FWE, within-basis component+tail correction, and conservative correction across all 32 components × two tails (64). Per-task means and network-template/atlas correspondence are descriptive. No new PPI or social modulation of connectivity is tested.\n\nAnatomy and social-specificity plots publish before whole-brain processing. They reuse existing participant task estimates: primary contrasts are UGR social−nonsocial fairness and Trust stranger−computer rec−def; friend comparisons/averages are secondary. All five original bars remain in source records. New contrasts are paired within participant; no cross-task coefficient difference is tested.\n\nAtlas target and independent Smith2009 RSN summaries exclude the defining ROI plus a predeclared 10-mm guard band to limit local spatial overlap. This guard band affects descriptive target ranking/template correspondence only, not the whole-brain inference mask. Labels are spatial correspondence, not forced network membership. Signed multiple-regression maps are partial associations, not causal effects.\n\nAll production participant/run/merged maps remain in the ignored output root. Small group maps, inference outputs, commands, masks, anatomy, summaries and failures publish to this record. All matching completed stage outputs are verified and reused on resume.\n')
        try:execute(args,record,root,state)
        except Exception as exc:state.update(state='wholebrain_dual_regression_failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        f.publish_stage(record,state,args.publish,'Publish functional whole-brain dual regression, TFCE networks and social task contrasts')
        print('Final state:',state['state'],'| record:',record,flush=True)
        raise SystemExit(0 if state['state']==SUCCESS else 1)


if __name__=='__main__':main()
