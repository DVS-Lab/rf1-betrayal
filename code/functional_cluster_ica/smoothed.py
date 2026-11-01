#!/usr/bin/env python3
"""Matched sm-5 functional-cluster arm; do not alter or rerun the sm-0 arm."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
import traceback
import uuid

import run as f
from rebuild_trust import recovery_fsf
import nibabel as nib
import numpy as np
from scipy import stats

ROOT=f.ROOT;RECORDS=f.RECORDS
POLICY='functional-sm5-matched-to-sm0-BOLD-design-mask-v1'
SUCCESS='matched_sm5_functional_cluster_complete_review_pending'
IGNORED={'fmri(outputdir)','fmri(smooth)','fmri(overwrite_yn)','feat_files(1)','confoundev_files(1)'}


def scientific_settings(settings):
    result={k:v for k,v in settings.items() if k not in IGNORED}
    for key in list(result):
        if re.fullmatch(r'fmri\(custom\d+\)',key):
            path=Path(result[key]);result[key]=f.g.core.sha(path) if path.is_file() else 'missing:'+str(path)
    return result


def source_signature(row,reference,source):
    feat=Path(row['active_feat']);design,ids=f.g.inputs.inspect_source(row,reference)
    frozen=f.read(source/'runs'/(row['key']+'.json'))['source_identities']
    f.g.core.require(ids==frozen,'Unsmooth fitted products changed since functional residual extraction')
    settings=f.g.core.fsl_settings(feat/'design.fsf');bold=f.n.feat_input(row,settings,reference)
    f.g.core.require(float(settings['fmri(smooth)'])==0,'Expected unsmoothed parent')
    events=[Path(v) for k,v in settings.items() if re.fullmatch(r'fmri\(custom\d+\)',k) and settings.get('fmri(shape'+re.search(r'\d+',k)[0]+')')=='3']
    f.g.core.require(all(p.is_file() for p in events),'Original seconds-based task EV missing')
    return dict(unsmoothed_files=ids,fsf=f.g.core.identity(feat/'design.fsf'),BOLD=f.g.core.identity(bold,False),
        events=[f.g.core.identity(p) for p in events],scientific_settings=scientific_settings(settings),
        nuisance_sha256=row['confounds_sha256'],TR=float(row['tr'])),settings,design


def check_matched(row,feat,signature,reference):
    """No candidate is reused on the basis of its name or smoothing label."""
    feat=Path(feat);settings=f.g.core.fsl_settings(feat/'design.fsf')
    f.g.core.require(float(settings['fmri(smooth)'])==5,'Candidate is not sm-5')
    f.g.core.require(scientific_settings(settings)==signature['scientific_settings'],'Candidate differs beyond spatial smoothing/output paths')
    bold=Path(settings['feat_files(1)'])
    f.g.core.require(f.g.core.identity(bold,False)==signature['BOLD'],'Candidate uses different BOLD')
    f.g.core.require(np.isclose(float(settings['fmri(tr)']),signature['TR'],rtol=0,atol=1e-5),'Candidate TR differs')
    old=Path(row['active_feat'])
    for name in ('design.mat','design.con'):
        a,b=f.g.core.fsl_matrix(old/name),f.g.core.fsl_matrix(feat/name)
        f.g.core.require(a.shape==b.shape and np.allclose(a,b,rtol=1e-8,atol=1e-8),'Matched full temporal model differs: '+name)
    a,b=f.g.core.load_numeric(row['confounds']),f.g.core.load_numeric(feat/'confoundevs.txt')
    f.g.core.require(np.array_equal(a,b),'Matched nuisance matrix differs')
    current=dict(row,active_feat=str(feat),original_feat=str(feat),source_fsf=str(feat/'design.fsf'),
        filtered_data=str(feat/'filtered_func_data.nii.gz'),native_residual=str(feat/'stats/res4d.nii.gz'),smoothing_mm=5)
    _,ids=f.g.inputs.inspect_source(current,reference)
    native=nib.load(current['native_residual'])
    f.g.core.require(f.g.core.same_grid(native,reference) and native.shape==(*reference.shape,row['nvol']) and np.isclose(native.header.get_zooms()[3],row['tr'],atol=1e-5),'Matched native residual shape/TR differs')
    return current,ids


def smoothed_fsf(text,row,target):
    settings=f.g.core.fsl_settings(Path(row['active_feat'])/'design.fsf')
    value=recovery_fsf(text,settings['feat_files(1)'],row['confounds'],target)
    value,count=re.subn(r'^set\s+fmri\(smooth\)\s+.*$','set fmri(smooth) 5',value,flags=re.M)
    f.g.core.require(count==1,'Missing/ambiguous smoothing setting')
    return value


def candidate_index(rows):
    index={r['key']:[] for r in rows}
    for row in rows:
        if row.get('manuscript_original_feat'):index[row['key']].append(Path(row['manuscript_original_feat']))
    # Exact previously recorded project paths only. No recursive disk search.
    parent=ROOT/'code/insula_ica/records/linux2-paired-group-ica-20261028T215022Z-679a3e/accepted-inputs.json'
    if parent.is_file():
        for row in f.read(parent):
            if row['key'] in index:index[row['key']].append(Path(row['active_feat']))
    for record in sorted((ROOT/'code/insula_ica/records').glob('linux2-ugr-TR-repair-*')):
        if (record/'state.json').is_file():
            state=f.read(record/'state.json');fit=state.get('fits',{}).get('corrected_sm5',{})
            if state.get('key') in index and fit.get('state')=='complete':index[state['key']].insert(0,Path(fit['feat']))
    return {k:list(dict.fromkeys(v)) for k,v in index.items()}


def fit_feat(row,root,record,source,reference,candidates):
    key=row['key'];signature,settings,design=source_signature(row,reference,source)
    marker=root/'FEAT-checkpoints'/(key+'.json')
    if marker.is_file():
        old=f.read(marker)
        f.g.core.require(old['signature']==signature,'Completed matched FEAT source changed')
        current,ids=check_matched(row,old['selected_feat'],signature,reference)
        f.g.core.require(ids==old['outputs'],'Completed matched FEAT products changed')
        f.save(record/'FEAT-snapshot/runs'/(key+'.json'),dict(source_identities=ids))
        return current,dict(reused_completed_feat=True,selected_feat=current['active_feat'],origin=old['origin'])
    rejected=[]
    for candidate in candidates:
        try:current,ids=check_matched(row,candidate,signature,reference)
        except Exception as exc:rejected.append(dict(path=str(candidate),reason=str(exc)));continue
        value=dict(signature=signature,selected_feat=current['active_feat'],outputs=ids,origin='verified_existing_sm5',rejected_candidates=rejected)
        f.save(marker,value);f.save(record/'FEAT-snapshot/runs'/(key+'.json'),dict(source_identities=ids))
        f.save(record/'FEAT-selection'/(key+'.json'),value)
        return current,dict(reused_completed_feat=True,selected_feat=current['active_feat'],origin=value['origin'])
    target=root/'feat'/(key+'_sm-5.feat');fsf=root/'fsf'/(key+'_sm-5.fsf')
    text=smoothed_fsf((Path(row['active_feat'])/'design.fsf').read_text(),row,target)
    if fsf.is_file():f.g.core.require(fsf.read_text()==text,'Prepared matched FSF changed')
    else:fsf.write_text(text)
    # Verify the generated temporal model BEFORE starting a new FEAT job.
    with (root/'logs'/(key+'.model.log')).open('w') as log:subprocess.run(['feat_model',str(fsf.with_suffix('')),row['confounds']],stdout=log,stderr=subprocess.STDOUT,check=True)
    for name,expected in (('mat',Path(row['active_feat'])/'design.mat'),('con',Path(row['active_feat'])/'design.con')):
        a,b=f.g.core.fsl_matrix(fsf.with_suffix('.'+name)),f.g.core.fsl_matrix(expected)
        f.g.core.require(a.shape==b.shape and np.allclose(a,b,rtol=1e-8,atol=1e-8),'Generated sm5 temporal design differs before FEAT')
    if target.exists():
        try:current,ids=check_matched(row,target,signature,reference)
        except Exception:
            processes=subprocess.check_output(['ps','-u',str(os.getuid()),'-o','args='],text=True)
            f.g.core.require(not any(str(fsf) in line and ('feat' in line or 'fsl' in line) for line in processes.splitlines()),'An earlier FEAT process still uses this FSF; wait before resuming')
            f.preserve(target);current=None
    else:current=None
    if current is None:
        with (root/'logs'/(key+'.feat.log')).open('w') as log:subprocess.run(['feat',str(fsf)],stdout=log,stderr=subprocess.STDOUT,check=True)
        current,ids=check_matched(row,target,signature,reference)
    f.g.core.require(source_signature(row,reference,source)[0]==signature,'Unsmoothed source changed during matched FEAT')
    value=dict(signature=signature,selected_feat=current['active_feat'],outputs=ids,origin='new_matched_sm5_FEAT',rejected_candidates=rejected)
    f.save(marker,value);f.save(record/'FEAT-snapshot/runs'/(key+'.json'),dict(source_identities=ids));f.save(record/'FEAT-selection'/(key+'.json'),value)
    return current,dict(reused_completed_feat=False,selected_feat=current['active_feat'],origin=value['origin'])


def load_task_cache(row,source_root,reference,task_source_output):
    model,cs,labels,_,anchors=f.task_context(row,reference,task_source_output)
    folder=source_root/'task_values'/row['key'];old=f.read(folder/'extraction.json');path=folder/'task_BOLD.npy'
    f.g.core.require(old['signature']==dict(anchors,mask_sha256=f.MASK_SHA) and f.g.core.identity(path)==old['output'],'Unsmoothed raw-task cache changed')
    values=np.load(path);f.verify_anchors(anchors)
    return values,model,cs,labels,anchors


def matched_plots(groups,means,existing,fits,folder):
    import matplotlib.pyplot as plt
    refs=[[r for r in means if r['series']=='functional_mean'],[r for r in existing if r['roi']=='left' and r['series']=='insula_mean']]
    for item in fits.values():
        if item['state']!='complete':continue
        mode=item['mode'];dim=item['dimensionality'];actual=item['actual_components']
        for start in range(1,actual+1,8):
            ks=list(range(start,min(start+8,actual+1)))
            blocks=refs+[[r for r in groups if r['mode']==mode and int(r['dimensionality'])==dim and int(r['component'])==k] for k in ks]
            fig,ax=plt.subplots(figsize=(max(12,len(blocks)*1.4),5.8))
            f.draw_bars(ax,blocks,['Functional mean','Left insula mean',*['IC'+str(k) for k in ks]],f'Matched 5 mm discovery | {mode} | requested {dim}, estimated {actual}')
            ax.legend(ncol=3,loc='upper center',bbox_to_anchor=(.5,-.12),frameon=False,fontsize=8)
            fig.text(.5,.012,'Same raw task-input BOLD, design, nuisance, subjects and mask as sm-0. Means reused. Exploratory task-selected ROI; pointwise 95% CIs.',ha='center',fontsize=8)
            fig.tight_layout(rect=(0,.06,1,1));fig.savefig(folder/f'sm5_{mode}_dim-{dim}_task_bars_page-{(start-1)//8+1}.png',dpi=150);plt.close(fig)


def compare_smoothing(source,state,mask,record,new_participants):
    old=f.read(source/'state.json');comparisons=[]
    for key,second in state['fits'].items():
        first=old['fits'].get(key,{})
        if second.get('state')!='complete' or first.get('state')!='complete':continue
        a,b=[nib.load(r['oriented_maps']).get_fdata()[mask] for r in (first,second)]
        comparisons.append(dict(mode=second['mode'],dimensionality=second['dimensionality'],**f.g.compare_maps(a,b,labels=('sm0','sm5'))))
    f.save(record/'smoothing_spatial_correspondence.json',comparisons)
    previous=source/'functional_IC_participant_effects.tsv'
    if not previous.is_file():return
    before=f.g.core.read_table(previous)
    def index(rows):return {(r['mode'],int(r['dimensionality']),int(r['component']),r['contrast'],r['subject']):float(r['effect']) for r in rows}
    a,b=index(before),index(new_participants);out=[]
    subjects=sorted({r['subject'] for r in before})
    for comparison in comparisons:
        mode=comparison['mode'];dim=comparison['dimensionality']
        for match in comparison['matches']:
            sign=match['sm5_sign_for_sm0_comparison']
            for label in f.n.BAR_LABELS:
                x=np.array([a[mode,dim,match['sm0_component'],label,s] for s in subjects]);y=np.array([sign*b[mode,dim,match['sm5_component'],label,s] for s in subjects])
                delta=y-x;se=float(stats.sem(delta));ci=float(stats.t.ppf(.975,len(subjects)-1)*se)
                out.append(dict(mode=mode,dimensionality=dim,sm0_component=match['sm0_component'],sm5_component=match['sm5_component'],spatial_r=match['spatial_r'],sm5_sign_for_comparison=sign,contrast=label,n=len(subjects),sm0_mean=float(x.mean()),aligned_sm5_mean=float(y.mean()),paired_sm5_minus_sm0=float(delta.mean()),ci_low=float(delta.mean()-ci),ci_high=float(delta.mean()+ci),inference='Descriptive pointwise paired-participant CI; spatial matching only, no significance selection'))
    f.g.core.write_table(record/'smoothing_paired_task_comparisons.tsv',out)


def execute(args,record,root,state):
    source=Path(state['unsmoothed_record']);old=f.read(source/'state.json');source_root=Path(old['output_root'])
    f.g.core.require(old['policy']==f.POLICY and len(old.get('residual_runs',{}))==490 and all(v['state']=='complete' for v in old['residual_runs'].values()),'Unsmooth functional preparation incomplete; do not rerun or substitute it')
    rows=f.read(source/'accepted-inputs.json')
    cohort=set((ROOT/'code/sublist_n132.txt').read_text().split())
    f.g.core.require(len(rows)==len({r['key'] for r in rows})==490 and len(cohort)==132 and {r['subject'] for r in rows}==cohort and Counter(r['task'] for r in rows)==dict(trust=241,ugr=249),'Matched cohort/run identities changed')
    reference=nib.load(Path(rows[0]['active_feat'])/'mean_func.nii.gz')
    mask_image=nib.load(source_root/'masks'/(f.ROI+'.nii.gz'));mask=mask_image.get_fdata()>0
    original=nib.load(source_root/'masks/original_functional_cluster.nii.gz').get_fdata()>0
    f.g.core.require(int(original.sum())==127 and np.all(mask<=original) and mask.sum()>=.8*127,'Invalid frozen sm0 functional mask')
    task_source=Path(old['source_record']);task_output=Path(f.read(task_source/'state.json')['output_root'])
    shared=dict(policy=POLICY,unsmoothed_record=str(source),source_registry=f.g.core.identity(source/'accepted-inputs.json'),original_mask_sha256=f.MASK_SHA,
        fixed_ICA_mask=f.g.core.identity(source_root/'masks'/(f.ROI+'.nii.gz')),source_means=f.g.core.identity(source/'functional_mean_group_effects.tsv'),
        code_sha256={Path(m.__file__).name:f.g.core.sha(m.__file__) for m in (f,sys.modules[__name__],f.n,f.d.pooled,f.g,f.g.inputs,f.retry)})
    f.d.freeze_json(record/'provenance.json',shared);f.d.freeze_json(record/'parent-accepted-inputs.json',rows)
    software=f.g.core.software()
    f.g.core.require(software['fsl']==f.read(source/'software.json')['fsl'],'FSL version differs from the unsmoothed arm')
    f.save(record/'software.json',software)
    f.d.freeze_json(record/'policy.json',dict(smoothing_mm=5,paired_with_smoothing_mm=0,cohort_n=132,accepted_runs=490,dimensions=list(f.DIMENSIONS),residual_modes=list(f.g.inputs.MODES),variance_normalization=True,temporal_concatenation=True,physically_masked_4D=True,same_frozen_ICA_mask=True,same_full_temporal_model=True,same_corrected_TRs=True,raw_task_BOLD_and_means_reused=True,no_upstream_preprocessing=True))
    for src,name in ((source/'functional_mean_group_effects.tsv','functional_mean_group_effects.tsv'),(source/'reused_insula_group_effects.tsv','reused_insula_group_effects.tsv')):f.d.copy_frozen(src,record/name)
    f.g.freeze_nifti(root/'masks'/(f.ROI+'.nii.gz'),mask.astype('uint8'),reference.affine)
    f.d.copy_frozen(source_root/'masks'/(f.ROI+'.nii.gz'),record/'primary_maps/ICA_common_functional_cluster.nii.gz')
    f.d.copy_frozen(source/'primary_maps/original_functional_cluster.nii.gz',record/'primary_maps/original_functional_cluster.nii.gz')
    gate=f.g.inputs.MemoryGate();lock=threading.Lock();candidates=candidate_index(rows)
    def update(section,key,value):
        with lock:
            state[section][key]=value;f.save(record/'state.json',state)
            print(f'{section}: {sum(v["state"]=="complete" for v in state[section].values())}/{8 if section=="fits" else 490} | {key}: {value["state"]}',flush=True)
    def worker(row,section,op):
        try:
            with gate.reserve(2):result,extra=op(row)
            status=dict(state='complete',**extra)
            error=record/'errors'/(section+'_'+row['key']+'.txt')
            if error.is_file():error.rename(record/'attempts'/(error.stem+'-'+uuid.uuid4().hex[:6]+'.txt'))
        except Exception as exc:
            result=None;status=dict(state='failed',error=str(exc));(record/'errors'/(section+'_'+row['key']+'.txt')).write_text(traceback.format_exc())
        update(section,row['key'],status);return row['key'],result
    def batch(section,op):
        pilots=[rows[0],next(r for r in rows if r['key']=='sub-11031_ses-01_task-ugr_run-1'),next(r for r in rows if r['key']=='sub-10817_ses-01_task-ugr_run-1')]
        results=dict(worker(r,section,op) for r in pilots);f.g.core.require(all(v is not None for v in results.values()),section+' pilot failed; no expansion')
        state['pilot_passed']=True;f.save(record/'state.json',state)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for future in as_completed([pool.submit(worker,r,section,op) for r in rows if r['key'] not in results]):
                key,value=future.result();results[key]=value
        f.g.core.require(all(v is not None for v in results.values()),section+' incomplete; completed work is preserved')
        return results
    state['state']='verifying_or_fitting_matched_sm5_FEAT';f.save(record/'state.json',state)
    selected=batch('feat_runs',lambda row:fit_feat(row,root,record,source,reference,candidates[row['key']]))
    sm5rows=[selected[r['key']] for r in rows];f.d.freeze_json(record/'accepted-inputs.json',sm5rows)
    f.publish_stage(record,state,args.publish,'Prepare matched 5 mm FEAT inputs for functional-cluster ICA')
    def residual(row):
        qc,reused=f.prepare_residuals(row,root,record,original,reference,record/'FEAT-snapshot')
        return qc,dict(residual_arrays_reused=reused)
    # The residual helper accepts these verified sm5 snapshots; it performs only
    # existing-PE algebra and native-residual extraction, without fitting FILM.
    saved_rows=rows;rows=sm5rows
    state['state']='matched_sm5_residual_extraction';f.save(record/'state.json',state)
    batch('residual_runs',residual)
    common=np.ones(127,bool)
    for row in rows:common &= np.load(root/'run_arrays'/(row['key']+'_valid.npy'))
    covered=bool(common[mask[original]].all())
    f.save(record/'coverage.json',dict(original_voxels=127,fixed_sm0_ICA_voxels=int(mask.sum()),sm5_common_valid_voxels=int(common.sum()),all_sm0_ICA_voxels_valid_in_sm5=covered,mask_identical_between_arms=True,no_erosion_or_exclusions=True))
    f.g.core.require(covered,'Matched sm5 data do not cover every frozen sm0 ICA voxel; no erosion or sm0 rerun performed')
    state['state']='building_physically_masked_sm5_blocks';f.save(record/'state.json',state)
    slices=f.build_sessions(rows,root,mask,original,reference);masks={f.ROI:mask}
    plans=[(m,f.ROI,d) for m in f.g.inputs.MODES for d in f.DIMENSIONS]
    f.save(record/'ica-plan.json',[dict(mode=m,roi=r,dimensionality=d,first_attempt_argv=f.g.melodic_command(m,r,d,root,output_folder=root/'ica_retries'/m/r/f'dim-{d}'/f'seed-{f.retry.SEEDS[0]}_maxit-{f.retry.MAXIT}',seed=f.retry.SEEDS[0],maxit=f.retry.MAXIT,variance_normalization=True,preserve_timing=True)) for m,r,d in plans])
    def fit(plan):
        mode,roi,dim=plan;key=f'{mode}/{roi}/dim-{dim}'
        try:
            if dim==0 and Path('/proc/meminfo').is_file():
                total=next(int(x.split()[1])/1024**2 for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemTotal:'))
                f.g.core.require(total>=40,'Automatic PCA needs 32 GiB plus 8 GiB headroom')
            with gate.reserve(32 if dim==0 else 2):value=f.retry.retry_fit(plan,root,masks,slices,reference,record,variance_normalization=True,preserve_timing=True,threads=8 if dim==0 else 4)
            f.d.copy_frozen(Path(value['oriented_maps']),record/'primary_maps'/f'sm5_{mode}_dim-{dim}.nii.gz')
            error=record/'errors'/(key.replace('/','_')+'.txt')
            if error.is_file():error.rename(record/'attempts'/(error.stem+'-'+uuid.uuid4().hex[:6]+'.txt'))
        except Exception as exc:value=dict(state='failed',mode=mode,roi=roi,dimensionality=dim,error=str(exc));(record/'errors'/(key.replace('/','_')+'.txt')).write_text(traceback.format_exc())
        update('fits',key,value)
    state['state']='six_fixed_matched_sm5_ICA_fits';f.save(record/'state.json',state)
    with ThreadPoolExecutor(max_workers=6) as pool:
        for future in as_completed([pool.submit(fit,p) for p in plans if p[2]!=0]):future.result()
    state['state']='two_serial_automatic_matched_sm5_ICA_fits';f.save(record/'state.json',state)
    for plan in plans:
        if plan[2]==0:fit(plan)
    bases={};map_ids={}
    for key,value in state['fits'].items():
        if value['state']!='complete':continue
        data=nib.load(value['oriented_maps']).get_fdata();f.g.core.require(np.count_nonzero(data[~mask])==0,'Nonzero sm5 component map outside fixed mask')
        bases[key]=(f.spatial_weights(data[mask]),value);map_ids[key]=f.g.core.identity(value['oriented_maps'])
        for start in range(0,value['actual_components'],12):f.g.component_figure(data[...,start:start+12],mask,reference.get_fdata(),reference.affine,record/'figures'/f'sm5_{value["mode"]}_dim-{value["dimensionality"]}_spatial_page-{start//12+1}.png',f'Matched 5 mm functional cluster | {key}',first_component=start+1)
    f.g.core.require(bases,'No sm5 ICA converged; completed FEAT/residuals retained')
    f.save(record/'task-map-identities.json',map_ids)
    def component(row):
        values,model,cs,labels,anchors=load_task_cache(row,source_root,reference,task_output);values=values[mask[original]].astype(float);effects=[];reused=[]
        for key,(weights,value) in bases.items():
            mode=value['mode'];dim=value['dimensionality'];layout=[dict(mode=mode,roi=f.ROI,dimensionality=dim,series='component',component=k) for k in range(1,value['actual_components']+1)]
            result,old=f.fit_response(row,root/'component_fits'/mode/f'dim-{dim}'/row['key'],dict(shared,anchors=anchors,maps={key:map_ids[key]},layout=layout),(weights@values).T,layout,model,cs,labels)
            effects.extend(result);reused.append(old)
        return effects,dict(raw_extraction_reused=True,basis_fits_reused=sum(reused),basis_fits_total=len(reused))
    rows=saved_rows;state['state']='matched_sm5_component_task_fits';f.save(record/'state.json',state)
    results=batch('component_runs',component);effects=[r for key in sorted(results) for r in results[key]]
    participant,groups=f.d.aggregate(effects,sorted({r['subject'] for r in rows}))
    for name,values in [('functional_IC_run_effects.tsv',effects),('functional_IC_participant_effects.tsv',participant),('functional_IC_group_effects.tsv',groups)]:f.g.core.write_table(record/name,values)
    means=f.g.core.read_table(record/'functional_mean_group_effects.tsv');existing=f.g.core.read_table(record/'reused_insula_group_effects.tsv')
    matched_plots(groups,means,existing,state['fits'],record/'figures');compare_smoothing(source,state,mask,record,participant)
    state.update(state=SUCCESS if len(state['fits'])==8 and all(v['state']=='complete' for v in state['fits'].values()) else 'matched_sm5_partial_ICA_review_errors',ICA_complete=len(bases),existing_unsmoothed_ICA_reruns=0,mean_fits_rerun=0)
    f.save(record/'state.json',state)


def latest_source():
    choices=sorted(p for p in RECORDS.glob('linux2-functional-cluster-*') if (p/'state.json').is_file())
    f.g.core.require(choices,'No existing unsmoothed functional-cluster run')
    return choices[-1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--unsmoothed-record',type=Path);p.add_argument('--resume-record',type=Path);p.add_argument('--jobs',type=int,default=64);p.add_argument('--publish',action='store_true');args=p.parse_args()
    f.g.core.require(1<=args.jobs<=64,'Expected 1–64 RAM-gated FEAT/preparation/task slots')
    source=args.unsmoothed_record or latest_source()
    if args.resume_record is None:
        for candidate in sorted(RECORDS.glob('linux2-functional-sm5-*'),reverse=True):
            if (candidate/'state.json').is_file():
                old=f.read(candidate/'state.json')
                if old.get('policy')==POLICY and Path(old['unsmoothed_record']).resolve()==source.resolve():args.resume_record=candidate;break
    output=ROOT/'code/insula_ica/output';output.mkdir(exist_ok=True)
    with (output/'functional-sm5.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:print('Matched sm5 submitter already active.');return
        if args.resume_record:
            record=args.resume_record;state=f.read(record/'state.json');root=Path(state['output_root']);source=Path(state['unsmoothed_record'])
            f.g.core.require(record.resolve().is_relative_to(RECORDS.resolve()) and root.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('functional-sm5-'),'Invalid sm5 resume root')
        else:
            stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6];record=RECORDS/('linux2-functional-sm5-'+stamp);root=ROOT/'derivatives/insula_ica'/('functional-sm5-'+stamp)
            record.mkdir();root.mkdir(parents=True);state=dict(policy=POLICY,state='waiting_for_unsmoothed_runner',unsmoothed_record=str(source),output_root=str(root),feat_runs={},residual_runs={},fits={},component_runs={},pilot_passed=False)
        for name in ('FEAT-selection','FEAT-snapshot/runs','runs','figures','primary_maps','attempts','errors'):(record/name).mkdir(parents=True,exist_ok=True)
        for name in ('FEAT-checkpoints','feat','fsf','logs','run_arrays','masks'):(root/name).mkdir(exist_ok=True)
        if (record/'failure.txt').is_file():(record/'failure.txt').rename(record/'attempts'/('failure-'+uuid.uuid4().hex[:6]+'.txt'))
        state.pop('error',None);state.update(state='waiting_for_unsmoothed_runner',jobs=args.jobs,runner_pid=os.getpid(),code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip());f.save(record/'state.json',state)
        print('Waiting for the existing unsmoothed runner to release its lock; it will not be interrupted.',flush=True)
        with (output/'functional-cluster.lock').open('a') as parent_lock:
            fcntl.flock(parent_lock,fcntl.LOCK_SH)
            # Keep a shared lock throughout this arm: it prevents a competing
            # sm0 resume from changing inputs, while allowing this reader.
            (record/'report.md').write_text('# Matched 5 mm functional-cluster arm\n\nSame source BOLD, corrected run-specific TR, full temporal model/contrasts, original nuisance matrix, cohort, fixed ICA mask and discovery settings as the sm0 arm. Only FEAT spatial smoothing is changed to 5 mm. Reuse exact verified sm5 fits; rebuild only missing/mismatched first-level fits in this separate root. No fMRIPrep/TEDANA or unsmoothed rerun. Smooth whole-brain BOLD in FEAT before masking; do not smooth an already masked ROI.\n\nBoth retained joint BOLD-scale and native FILM residuals; physical 4D masking; 3/5/7/automatic concat ICA with native variance normalization, identical balancing and serialized automatic PCA. Require sm5 support at every frozen sm0 ICA voxel; no erosion, voxel/subject selection or mask substitution. All masks/inputs/settings and candidate rejection reasons are recorded.\n\nRaw FEAT-input task BOLD, mean estimates and anatomical-insula results are reused. New sm5 component maps undergo the same spatial regression on that identical task-preserved input and pooled-UGR/full task+nuisance standardized OLS. Five bars per component, pointwise participant CIs. Spatial-only cross-smoothing matching and paired task differences are descriptive; components are not assumed to correspond by index. This remains exploratory task-selected ROI characterization.\n')
            try:execute(args,record,root,state)
            except Exception as exc:state.update(state='matched_sm5_failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        f.publish_stage(record,state,args.publish,'Publish matched smoothed functional-cluster ICA and paired smoothing comparisons')
        print('Final state:',state['state'],'| record:',record,flush=True)
        raise SystemExit(0 if state['state']==SUCCESS else 1)


if __name__=='__main__':main()
