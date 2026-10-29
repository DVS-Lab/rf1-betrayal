#!/usr/bin/env python3
"""Matched confounds-only/full-model OLS sensitivity on retained sm-5 data."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import fcntl
import json
from pathlib import Path
import subprocess
import threading
import traceback
import uuid

import nibabel as nib
import numpy as np

import run_group_ica as g
import resume_group_ica as resume
from run_unsmoothed_feat import SOURCE


MODES=('confounds_only_ols','full_model_ols')


def ols_pair(values,design,task_columns):
    y=np.asarray(values,dtype=float);x=np.asarray(design,dtype=float)
    g.core.require(y.ndim==2 and x.ndim==2 and y.shape[1]==x.shape[0] and 0<task_columns<x.shape[1],'OLS sensitivity dimensions disagree')
    g.core.require(np.isfinite(y).all() and np.isfinite(x).all(),'Nonfinite OLS sensitivity input')
    y=y-y.mean(axis=1,keepdims=True);out={}
    for mode,model in zip(MODES,(x[:,task_columns:],x)):
        q=g.core.basis(model);out[mode]=y-(y@q)@q.T
        g.core.require(np.linalg.norm(out[mode]@q)<=1e-8*max(np.linalg.norm(y),1),'Nuisance/full-model OLS orthogonality failed')
    return out


def prepare(row,root,record,union,reference):
    feat=Path(row['active_feat']);filtered=Path(row['filtered_data']);designpath=feat/'design.mat'
    settings=g.core.fsl_settings(feat/'design.fsf' if (feat/'design.fsf').is_file() else Path(row['source_fsf']))
    g.core.require(float(settings['fmri(smooth)'])==5 and settings['fmri(temphp_yn)']=='0','Expected original sm-5/no-high-pass source')
    g.core.require(g.core.sha(row['confounds'])==row['confounds_sha256'],'Saved TEDANA/nuisance matrix changed')
    design=g.core.fsl_matrix(designpath);nuisance=g.core.load_numeric(row['confounds']);taskcols=int(settings['fmri(evs_real)'])
    g.core.require(design.shape==(row['nvol'],taskcols+nuisance.shape[1]),'Cannot identify the exact appended nuisance columns')
    signature=dict(filtered_data=g.core.identity(filtered,False),design=g.core.identity(designpath),nuisance_sha256=row['confounds_sha256'],task_columns=taskcols)
    folder=root/'run_arrays';marker=folder/(row['key']+'_complete.json')
    if marker.is_file():
        old=json.loads(marker.read_text());g.core.require(old['signature']==signature,'Task-preserved source changed')
        g.core.require(all(g.core.identity(p['path'])==p for p in old['arrays']),'Task-preserved arrays changed')
        g.save(record/'runs'/(row['key']+'.json'),old['qc']);return old['qc']
    image=nib.load(filtered);g.core.require(g.core.same_grid(image,reference) and image.shape==(*reference.shape,row['nvol']),'Task-preserved input geometry mismatch')
    y=np.asarray(image.dataobj,dtype='float32')[union].astype(float)
    arrays=ols_pair(y,design,taskcols)
    for values in arrays.values():g.core.require((values.std(axis=1)>1e-8).all(),'Constant sensitivity residual inside fixed shared mask')
    paths=[folder/(row['key']+'_'+mode+'.npy') for mode in MODES]
    g.core.require(not any(p.exists() for p in paths),'Partial task-preserved arrays exist; preserve them')
    for mode,path in zip(MODES,paths):np.save(path,arrays[mode].astype('float32'))
    nu=arrays[MODES[0]];full=arrays[MODES[1]]
    qc=dict(key=row['key'],task_columns=taskcols,nuisance_columns=nuisance.shape[1],
        nuisance_rank=g.core.basis(design[:,taskcols:]).shape[1],full_model_rank=g.core.basis(design).shape[1],
        conditional_task_fraction_of_nuisance_only_residual_energy=float(np.sum((nu-full)**2)/np.sum(nu**2)),
        exact_fitted_nuisance_columns=True,task_regressors_removed_from_confounds_only=False,source_signature=signature)
    g.core.require(g.core.identity(filtered,False)==signature['filtered_data'] and g.core.identity(designpath)==signature['design'],'Sensitivity source changed during preparation')
    g.save(record/'runs'/(row['key']+'.json'),qc);g.save(marker,dict(signature=signature,qc=qc,arrays=[g.core.identity(p) for p in paths]))
    return qc


def run(args,root,record,state):
    source_state=json.loads((args.source_record/'state.json').read_text());parent=Path(source_state['output_root'])
    _,prep_source,_,masks,slices,reference=resume.verify_prepared(args.source_record,parent,record)
    rows=json.loads((prep_source/'accepted-inputs.json').read_text())
    g.save(record/'policy.json',dict(cohort_n=132,trust_runs=241,ugr_runs=249,modes=list(MODES),smoothing_mm=5,
        components=5,rois=list(g.core.ROIS),fitting='Matched OLS, existing complete design versus its exact appended nuisance block.',
        nuisance_only_interpretation='Retains task variance outside the nuisance span; task-correlated nuisance variance is removed.',
        native_melodic_variance_normalization=True,global_input_amplitude_multiplier=1000,
        no_FEAT_or_fMRIPrep=True,no_new_temporal_filter=True,no_additional_smoothing=True))
    g.save(record/'software.json',g.core.software())
    (root/'masks').mkdir(exist_ok=True)
    for roi,mask in masks.items():g.freeze_nifti(root/'masks'/(roi+'.nii.gz'),mask.astype('uint8'),reference.affine)
    state.update(state='preparing_matched_OLS_sensitivity',preparation_source_record=str(prep_source));g.save(record/'state.json',state)
    gate=g.inputs.MemoryGate();lock=threading.Lock()
    def worker(row):
        try:
            with gate.reserve(1):prepare(row,root,record,masks['bilateral'],reference)
            result=dict(state='complete')
        except Exception as exc:
            result=dict(state='failed',error=str(exc));(record/'errors'/(row['key']+'.txt')).write_text(traceback.format_exc())
        with lock:
            state['runs'][row['key']]=result;g.save(record/'state.json',state)
            if len(state['runs'])%20==0 or result['state']=='failed':print(f'Task-preserved inputs {len(state["runs"])}/490 | {row["key"]}: {result["state"]}',flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(worker,row) for row in rows]):future.result()
    g.core.require(all(v['state']=='complete' for v in state['runs'].values()),'Task-preserved input failures; no sensitivity ICA started')
    state['state']='building_task_preserved_masked_4d_inputs';g.save(record/'state.json',state)
    slices,_=g.sessions(rows,root,record,masks,masks['bilateral'],reference,modes=MODES,amplitude=1000)
    state['state']='running_6_task_preserved_sensitivity_fits';g.save(record/'state.json',state)
    def fit(mode,roi):
        key=f'{mode}/{roi}/dim-5'
        try:
            with gate.reserve(2):result=resume.retry_fit((mode,roi,5),root,masks,slices,reference,record,variance_normalization=True)
        except Exception as exc:
            result=dict(state='failed',error=str(exc));(record/'errors'/(key.replace('/','_')+'.txt')).write_text(traceback.format_exc())
        with lock:state['fits'][key]=result;g.save(record/'state.json',state)
        print('Task-preserved ICA '+key+': '+result['state'],flush=True)
    with ThreadPoolExecutor(max_workers=6) as pool:
        for future in as_completed([pool.submit(fit,mode,roi) for mode in MODES for roi in g.core.ROIS]):future.result()
    g.core.require(all(v['state']=='complete' for v in state['fits'].values()),'Some task-preserved sensitivity fits remain incomplete')
    comparisons=[]
    for roi in g.core.ROIS:
        a,b=[state['fits'][f'{mode}/{roi}/dim-5'] for mode in MODES]
        data=[nib.load(r['oriented_maps']).get_fdata() for r in (a,b)]
        comparisons.append(dict(roi=roi,first='confounds_only_ols',second='full_model_ols',**g.compare_maps(data[0][masks[roi]],data[1][masks[roi]],labels=MODES)))
        for mode,maps in zip(MODES,data):g.component_figure(maps,masks[roi],reference.get_fdata(),reference.affine,
            record/'figures'/f'{mode}_{roi}_dim-5.png',f'{mode} | {roi} | 5 components | native MELODIC VN')
    g.save(record/'task_removal_spatial_comparisons.json',comparisons)
    state['state']='task_preserved_sensitivity_complete_review_pending';g.save(record/'state.json',state)
    (record/'README.md').write_text('# Task-preserved ICA sensitivity\n\n'
        'N=132, Trust241/UGR249; existing sm-5 filtered data and exact saved full-model/nuisance designs. '
        'Matched OLS: nuisance-only versus full task+nuisance. Five components for each left/right/bilateral ROI; six fits. '
        'Native MELODIC variance normalization enabled in both arms; physically masked balanced inputs uniformly multiplied by 1000. '
        'This holds the fitting convention fixed within the task-removal contrast, while the original GLS-based series remains unchanged. '
        'Task-correlated nuisance variance is removed in both arms; the confounds-only arm retains remaining task variance. '
        'No fMRIPrep, FEAT, FILM, extra temporal filtering or spatial smoothing was performed. See figures and task_removal_spatial_comparisons.json.\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-record',type=Path,default=SOURCE);parser.add_argument('--resume-record',type=Path)
    parser.add_argument('--jobs',type=int,default=64);parser.add_argument('--publish',action='store_true');args=parser.parse_args()
    g.core.require(1<=args.jobs<=64,'Expected 1–64 preparation slots')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    if args.resume_record:
        previous=json.loads((args.resume_record/'state.json').read_text());root=Path(previous['output_root']).resolve();args.source_record=Path(previous['source_record'])
        g.core.require(root.is_relative_to((g.core.ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('task-preserved-ols-'),'Invalid sensitivity resume root')
    else:root=g.core.ROOT/'derivatives/insula_ica'/('task-preserved-ols-'+stamp)
    record=g.core.ROOT/'code/insula_ica/records'/('linux2-task-preserved-ica-'+stamp)
    lockpath=g.core.ROOT/'code/insula_ica/output/task-preserved-ica.lock';lockpath.parent.mkdir(parents=True,exist_ok=True)
    with lockpath.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        (root/'run_arrays').mkdir(parents=True,exist_ok=True)
        for name in ('runs','figures','errors','attempts'):(record/name).mkdir(parents=True,exist_ok=True)
        state=dict(state='preflight',runs={},fits={},output_root=str(root),source_record=str(args.source_record),
            code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=g.core.ROOT,text=True).strip())
        try:run(args,root,record,state)
        except Exception as exc:
            state.update(state='failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        g.save(record/'state.json',state)
        if args.publish:
            try:g.publish(record,'Compare task-preserved confounds-only and matched full-model ICA residuals')
            except Exception as exc:(record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
        print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state']=='task_preserved_sensitivity_complete_review_pending' else 1)


if __name__=='__main__':main()
