#!/usr/bin/env python3
"""Separate sm-0 first-level FEAT stream using the recorded BOLD/EV/confound inputs."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import threading
import traceback
import uuid

import nibabel as nib
import numpy as np

import run_group_ica as g
from rebuild_trust import recovery_fsf


SOURCE=g.core.ROOT/'code/insula_ica/records/linux2-paired-group-ica-20261028T215022Z-679a3e'


def unsmoothed_fsf(text,bold,confounds,target):
    updated=recovery_fsf(text,bold,confounds,target)
    updated,n=re.subn(r'^set\s+fmri\(smooth\)\s+.*$', 'set fmri(smooth) 0',updated,flags=re.M)
    g.core.require(n==1,'Expected exactly one saved smoothing setting')
    return updated


def plan(row,root,bold_override=None):
    feat=Path(row['active_feat']);fsf=feat/'design.fsf'
    if not fsf.is_file():fsf=Path(row['source_fsf'])
    settings=g.core.fsl_settings(fsf)
    g.core.require(float(settings['fmri(smooth)'])==5,'Expected the accepted sm-5 source setup')
    design=g.core.fsl_matrix(feat/'design.mat');contrasts=g.core.fsl_matrix(feat/'design.con')
    g.residuals.validate_model(settings,design,contrasts,row['task'])
    g.core.require(g.core.sha(row['confounds'])==row['confounds_sha256'],'Historical nuisance file changed')
    nuisance=g.core.load_numeric(row['confounds']);fitted=g.core.load_numeric(feat/'confoundevs.txt')
    g.core.require(nuisance.shape==fitted.shape and np.array_equal(nuisance,fitted),'Original nuisance columns changed')
    historical_bold=Path(settings['feat_files(1)'])
    if bold_override:
        g.core.require(row['task']=='ugr','Regenerated BOLD override is restricted to UGR')
        bold=Path(bold_override['path'])
        g.core.require(g.core.identity(bold,content=False)==bold_override['identity'],'Regenerated UGR BOLD changed')
    else:bold=historical_bold
    g.core.require(bold.is_absolute(),'Saved BOLD path is not absolute')
    events=[Path(settings[f'fmri(custom{k})']) for k in range(1,int(settings['fmri(evs_orig)'])+1)
        if settings.get(f'fmri(shape{k})')=='3']
    missing=[str(p) for p in (bold,*events) if not p.is_file()]
    target=root/'feat'/(row['key']+'_sm-0.feat')
    text=unsmoothed_fsf(fsf.read_text(),bold,row['confounds'],target)
    result=dict(row,source_sm5_feat=str(feat),source_settings=str(fsf),source_fsf_sha256=g.core.sha(fsf),
        source_bold=str(bold),source_bold_present=bold.is_file(),missing_inputs=missing,
        status='missing_source_inputs' if missing else 'ready',target=str(target),source_smoothing_mm=5,smoothing_mm=0,
        unsmoothed_fsf=str(root/'fsf'/(row['key']+('_sm-0_regenerated-ugr.fsf' if bold_override else '_sm-0.fsf'))))
    if bold_override:result.update(historical_bold=str(historical_bold),regenerated_bold_provenance=bold_override)
    if not missing:
        image=nib.load(bold);reference=nib.load(feat/'mean_func.nii.gz')
        g.core.require(g.core.same_grid(image,reference) and image.shape==(*reference.shape,row['nvol']),'Saved FEAT BOLD input has wrong grid/time count')
        result['source_bold_identity']=g.core.identity(bold,content=False)
        result['event_identities']=[g.core.identity(p) for p in events]
    path=Path(result['unsmoothed_fsf'])
    if path.is_file():g.core.require(path.read_text()==text,'Previously generated sm-0 FSF changed')
    else:path.write_text(text)
    return result


def fit(item,root,record):
    target=Path(item['target']);signature=dict(bold=item['source_bold_identity'],events=item['event_identities'],
        original_fsf_sha256=item['source_fsf_sha256'],confounds_sha256=item['confounds_sha256'],
        unsmoothed_fsf_sha256=g.core.sha(item['unsmoothed_fsf']))
    marker=target/'sm0_complete.json'
    if marker.is_file():
        old=json.loads(marker.read_text());g.core.require(old['signature']==signature,'Completed sm-0 inputs changed')
        g.core.require(all(g.core.identity(p['path'],content=p['sha256'] is not None)==p for p in old['outputs']),'Completed sm-0 outputs changed')
        g.save(record/'qc'/(item['key']+'.json'),old['qc']);return dict(state='complete',reused_completed_feat=True)
    g.core.require(not target.exists(),'Partial FEAT exists; preserve it for inspection: '+str(target))
    log=root/'logs'/(item['key']+'.txt')
    with log.open('w') as stream:subprocess.run(['feat',item['unsmoothed_fsf']],stdout=stream,stderr=subprocess.STDOUT,check=True)
    # Changing spatial smoothing must not alter the complete temporal design or nuisances.
    original=Path(item['source_sm5_feat'])
    for name in ('design.mat','design.con'):
        a,b=g.core.fsl_matrix(original/name),g.core.fsl_matrix(target/name)
        g.core.require(a.shape==b.shape and np.allclose(a,b,rtol=1e-8,atol=1e-8),'Unsmoothed complete model changed: '+name)
    a,b=g.core.load_numeric(item['confounds']),g.core.load_numeric(target/'confoundevs.txt')
    g.core.require(a.shape==b.shape and np.array_equal(a,b),'Unsmoothed nuisance file differs')
    settings=g.core.fsl_settings(target/'design.fsf');g.core.require(float(settings['fmri(smooth)'])==0,'FEAT smoothing is not zero')
    reference=nib.load(original/'mean_func.nii.gz');residual=nib.load(target/'stats/res4d.nii.gz')
    g.core.require(g.core.same_grid(reference,residual) and residual.shape==(*reference.shape,item['nvol']),'Unsmoothed residual geometry mismatch')
    mask=nib.load(target/'mask.nii.gz').get_fdata()>0
    values=np.asarray(residual.dataobj,dtype='float32');g.core.require(np.isfinite(values[mask]).all(),'Nonfinite unsmoothed residuals')
    valid=mask.copy();valid[mask]&=values[mask].std(axis=1)>1e-8
    atlas=Path(os.environ['FSLDIR'])/'data/atlases'
    masks,_,_=g.core.anatomical_masks(atlas/'HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz',atlas/'HarvardOxford-Cortical.xml',reference)
    dof=float((target/'stats/dof').read_text());g.core.require(dof>=37,'Insufficient sm-0 residual degrees of freedom')
    paths=[target/p for p in ('filtered_func_data.nii.gz','stats/res4d.nii.gz','stats/threshac1.nii.gz','stats/sigmasquareds.nii.gz',
        'design.mat','design.con','confoundevs.txt','mask.nii.gz','mean_func.nii.gz')]
    paths.extend(target/f'stats/pe{k+1}.nii.gz' for k in range(g.core.fsl_matrix(target/'design.mat').shape[1]))
    g.core.require(all(p.is_file() for p in paths),'Required sm-0 residual/PE/intermediate output missing')
    qc=dict(key=item['key'],smoothing_mm=0,complete_design_preserved=True,exact_original_nuisance_matrix=True,
        residual_dof=dof,finite_residuals=True,coverage={roi:float((m&valid).sum()/m.sum()) for roi,m in masks.items()},
        residual=str(target/'stats/res4d.nii.gz'),filtered_data=str(target/'filtered_func_data.nii.gz'),
        source_bold=item['source_bold'],source_bold_identity=item['source_bold_identity'],
        residual_scale='Native FILM voxelwise-prewhitened; no additional nuisance regression.',qc_review_pending=True)
    g.core.require(g.core.identity(item['source_bold'],content=False)==item['source_bold_identity'],'FEAT source BOLD changed during the fit')
    outputs=[g.core.identity(p,content=not p.name.endswith('.nii.gz')) for p in paths]
    g.save(record/'qc'/(item['key']+'.json'),qc);g.save(marker,dict(signature=signature,outputs=outputs,qc=qc))
    return dict(state='complete',reused_completed_feat=False)


def run(args,root,record,state):
    rows=json.loads((args.source_record/'accepted-inputs.json').read_text());g.core.require(len(rows)==490,'Expected the frozen 490 accepted runs')
    override_path=getattr(args,'bold_overrides',None)
    overrides=json.loads(override_path.read_text()) if override_path else {}
    g.core.require(set(overrides)<={r['key'] for r in rows if r['task']=='ugr'},'Unknown/non-UGR BOLD overrides')
    if override_path:
        state['bold_overrides']=str(override_path.resolve())
        g.save(record/'regenerated-ugr-bold-overrides.json',overrides)
    g.save(record/'software.json',g.core.software());plans=[]
    for row in rows:
        try:plans.append(plan(row,root,overrides.get(row['key'])))
        except Exception as exc:plans.append(dict(key=row['key'],task=row['task'],status='blocked',error=str(exc)))
    g.save(record/'input-plan.json',plans)
    ready=[p for p in plans if p['status']=='ready']
    state.update(state='sm0_input_audit_complete',input_status_counts=dict(Counter(p['status'] for p in plans)),jobs=args.jobs)
    g.save(record/'state.json',state)
    print('sm-0 exact-path input audit:',state['input_status_counts'],flush=True)
    for task in ('trust','ugr'):print(task,dict(Counter(p['status'] for p in plans if p['task']==task)),flush=True)
    gate=g.inputs.MemoryGate();lock=threading.Lock()
    if args.execute:
        state['state']='running_sm0_FEAT';g.save(record/'state.json',state)
        def worker(item):
            try:
                with gate.reserve(1):result=fit(item,root,record)
            except Exception as exc:
                result=dict(state='failed',error=str(exc));(record/'errors'/(item['key']+'.txt')).write_text(traceback.format_exc())
            with lock:
                state['runs'][item['key']]=result;g.save(record/'state.json',state)
                if len(state['runs'])%10==0 or result['state']=='failed':print(f'sm-0 FEAT {len(state["runs"])}/{len(ready)} | {item["key"]}: {result["state"]}',flush=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for future in as_completed([pool.submit(worker,p) for p in ready]):future.result()
        failed=sum(r['state']=='failed' for r in state['runs'].values())
        state['state']='sm0_FEAT_complete_QC_pending' if len(ready)==490 and failed==0 else 'sm0_partial_inputs_or_failed_fits'
    g.save(record/'state.json',state)
    (record/'README.md').write_text('# Separate unsmoothed FEAT stream\n\n'
        'Recorded task EVs and saved TEDANA+nuisance matrix; FEAT smoothing changes from 5 mm to 0 mm, with new output locations. '
        'Original BOLD is used unless an explicit, identity-checked regenerated UGR override is listed in regenerated-ugr-bold-overrides.json. '
        'No original wrapper, fMRIPrep, TEDANA or Warpkit was executed. No existing sm-5 fit was overwritten. '
        'All sm-0 filtered data, native residuals, autocorrelations and PEs are retained. '
        f'Worker queue slots: {args.jobs}, with 1 GiB reservations and 8 GiB headroom. '
        'See input-plan.json for exact saved BOLD paths and missing inputs. Available inputs are processed independently; '
        'missing pre-smoothing BOLD is never replaced with already-smoothed filtered_func_data. '
        'This record prepares the separate FEAT stream; registration/coverage and ICA review remain separate.\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-record',type=Path,default=SOURCE)
    parser.add_argument('--jobs',type=int,default=64)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--resume-record',type=Path)
    parser.add_argument('--bold-overrides',type=Path,help='Explicit identity/provenance registry for regenerated UGR BOLD only')
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args();g.core.require(1<=args.jobs<=64,'Expected 1–64 workers')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    if args.resume_record:
        previous=json.loads((args.resume_record/'state.json').read_text());root=Path(previous['output_root']).resolve()
        g.core.require(root.is_relative_to((g.core.ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('sm-0-feat-'),'Invalid sm-0 resume root')
        args.source_record=Path(previous['source_record'])
        if args.bold_overrides is None and previous.get('bold_overrides'):args.bold_overrides=Path(previous['bold_overrides'])
    else:root=g.core.ROOT/'derivatives/insula_ica'/('sm-0-feat-'+stamp)
    record=g.core.ROOT/'code/insula_ica/records'/('linux2-sm-0-feat-'+stamp)
    lockpath=g.core.ROOT/'code/insula_ica/output/sm-0-feat.lock';lockpath.parent.mkdir(parents=True,exist_ok=True)
    with lockpath.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for name in ('fsf','logs','feat'):(root/name).mkdir(parents=True,exist_ok=True)
        for name in ('qc','errors'):(record/name).mkdir(parents=True,exist_ok=True)
        state=dict(state='preflight',runs={},output_root=str(root),source_record=str(args.source_record),
            code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=g.core.ROOT,text=True).strip())
        try:run(args,root,record,state)
        except Exception as exc:
            state.update(state='failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        g.save(record/'state.json',state)
        if args.publish:
            try:g.publish(record,'Prepare separate unsmoothed sm-0 FEAT stream from recorded inputs')
            except Exception as exc:
                (record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
        print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state'] in ('sm0_FEAT_complete_QC_pending','sm0_input_audit_complete') else 1)


if __name__=='__main__':main()
