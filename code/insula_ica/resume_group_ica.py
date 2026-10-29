#!/usr/bin/env python3
"""Resume only unfinished ICA fits from verified, already prepared 4D inputs."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import re
import threading
import traceback
import uuid

import nibabel as nib
import numpy as np

import run_group_ica as g


SEEDS=tuple(g.SEED+k for k in range(5))
MAXIT=10000


def read(path):return json.loads(path.read_text())


def verify_prepared(source,root,record):
    previous=read(source/'state.json')
    prep_source=Path(previous.get('preparation_source_record',str(source)))
    prepared=read(prep_source/'state.json')
    g.core.require(len(prepared['runs'])==490 and all(r['state']=='complete' for r in prepared['runs'].values()),'All 490 runs must already be prepared')
    g.core.require(g.core.sha(root/'policy.json')==g.core.sha(prep_source/'policy.json'),'Prepared analysis policy changed')
    for key in prepared['runs']:
        checkpoint=read(root/'run_arrays'/(key+'_complete.json'))
        g.core.require(checkpoint['qc']==read(prep_source/'runs'/(key+'.json')),'Prepared run QC changed: '+key)
        g.core.require(all(g.core.identity(a['path'])==a for a in checkpoint['arrays']),'Prepared residual array changed: '+key)
    expected={(v['mode'],v['roi']):v['verified_inputs_sha256'] for v in read(prep_source/'physical-4d-masking.json')}
    verified={(mode,roi):g.verify_ica_inputs(mode,roi,root) for mode in g.inputs.MODES for roi in g.core.ROIS}
    g.core.require(verified==expected,'Previously verified physical 4D input collection changed')
    rows=read(prep_source/'accepted-inputs.json')
    reference=nib.load(Path(rows[0]['active_feat'])/'mean_func.nii.gz')
    geometry=read(prep_source/'crop.json')
    g.core.require(list(reference.shape)==geometry['full_shape'] and np.allclose(reference.affine,geometry['full_affine']),'Original ICA reference grid changed')
    slices=tuple(slice(a,b) for a,b in zip(geometry['crop_start'],geometry['crop_stop']))
    masks={roi:nib.load(root/'masks'/(roi+'.nii.gz')).get_fdata()>0 for roi in g.core.ROIS}
    for roi in g.core.ROIS:
        crop=nib.load(root/'masks'/(roi+'_crop.nii.gz'))
        full=nib.load(root/'masks'/(roi+'.nii.gz'))
        g.core.require(g.core.same_grid(full,reference) and np.array_equal(masks[roi][slices],crop.get_fdata()>0),'Full and cropped masks differ')
        translation=np.eye(4);translation[:3,3]=geometry['crop_start']
        g.core.require(np.allclose(reference.affine@translation,crop.affine),'ROI crop coordinates changed')
    g.save(record/'reuse-validation.json',dict(prepared_runs=490,verified_4d_inputs=1584,verified_input_collections={f'{m}/{r}':v for (m,r),v in verified.items()},
        preparation_source_record=str(prep_source),source_policy_sha256=g.core.sha(prep_source/'policy.json'),
        no_preprocessing=True,no_film=True,no_residual_reconstruction=True,no_4d_rewriting=True))
    return prepared,prep_source,verified,masks,slices,reference


def validate_completed(result,root,verified):
    mode,roi,dim=result['mode'],result['roi'],result['dimensionality']
    folder=Path(result['oriented_maps']).parent
    g.core.require(folder.resolve().is_relative_to(root.resolve()),'Completed fit lies outside this analysis root')
    saved=read(folder/'complete.json')
    default=root/'ica'/mode/roi/f'dim-{dim}'
    kwargs={} if folder==default else dict(output_folder=folder,seed=saved['seed'],maxit=saved['max_iterations'])
    g.core.require(saved['argv']==g.melodic_command(mode,roi,dim,root,**kwargs),'Completed ICA settings changed')
    g.core.require(saved['verified_inputs_sha256']==verified[mode,roi] and
        saved['mask_sha256']==g.core.sha(root/'masks'/(roi+'.nii.gz')),'Completed ICA inputs/mask changed')
    g.core.require(g.core.identity(saved['oriented_maps'])==saved['map_identity'],'Completed ICA map changed')
    return dict(saved,state='complete',reused_completed_ica=True)


def is_nonconvergence(folder):
    paths=[folder.parent/(folder.name+'.log'),folder/'log.txt']
    return not (folder/'melodic_IC.nii.gz').exists() and any(
        p.is_file() and re.search(r'No convergence',p.read_text(errors='replace'),re.I) for p in paths)


def record_attempt(record,folder,attempt):
    from collect_group_ica_diagnostics import tail
    prefix=('native_vn_' if attempt['variance_normalization'] else '')+f'{attempt["mode"]}_{attempt["roi"]}_dim-{attempt["dimensionality"]}_seed-{attempt["seed"]}'
    logs=[]
    for kind,path in (('stdout',folder.parent/(folder.name+'.log')),('melodic',folder/'log.txt')):
        if path.is_file():
            text,identity=tail(path);saved=record/'attempts'/(prefix+'_'+kind+'.txt');saved.write_text(text)
            logs.append(dict(**identity,saved=saved.name))
            converged=re.search(r'Convergence after\s+(\d+)\s+steps',text)
            if converged:attempt['converged_after_iterations']=int(converged.group(1))
    g.save(record/'attempts'/(prefix+'.json'),dict(attempt,logs=logs))


def retry_fit(plan,root,masks,slices,reference,record,on_attempt=None,*,variance_normalization=False):
    mode,roi,dim=plan
    for seed in SEEDS:
        folder=root/'ica_retries'/mode/roi/f'dim-{dim}'/f'seed-{seed}_maxit-{MAXIT}'
        note=folder.with_name(folder.name+'.attempt.json')
        kwargs=dict(output_folder=folder,seed=seed,maxit=MAXIT,variance_normalization=variance_normalization)
        if folder.exists() and not (folder/'complete.json').is_file():
            g.core.require(is_nonconvergence(folder),'Unverified partial retry exists; preserve it and inspect: '+str(folder))
            continue
        attempt=dict(mode=mode,roi=roi,dimensionality=dim,seed=seed,max_iterations=MAXIT,epsilon=0.00005,
            nonlinearity='pow3',variance_normalization=variance_normalization,argv=g.melodic_command(mode,roi,dim,root,**kwargs),state='running')
        folder.parent.mkdir(parents=True,exist_ok=True);g.save(note,attempt)
        if on_attempt:on_attempt(attempt)
        try:
            result=g.fit_one(mode,roi,dim,root,masks,slices,reference,**kwargs)
            attempt['state']='complete';g.save(note,attempt)
            record_attempt(record,folder,attempt)
            return dict(result,state='complete',reused_completed_ica=result.get('reused_completed_ica',False))
        except Exception as exc:
            nonconverged=is_nonconvergence(folder)
            attempt.update(state='nonconvergence' if nonconverged else 'failed',error=str(exc))
            g.save(note,attempt);record_attempt(record,folder,attempt)
            if not nonconverged:raise
            print(f'ICA {mode}/{roi}/dim-{dim}: seed {seed} did not converge',flush=True)
    raise ValueError('No convergence after five deterministic starts at 10000 iterations; inputs, dimension, pow3 objective and epsilon unchanged. All attempts preserved.')


def run(source,root,record,state,jobs):
    prepared,prep_source,verified,masks,slices,reference=verify_prepared(source,root,record)
    state.update(preparation_source_record=str(prep_source),runs=prepared['runs'],state='validating_completed_ICA_fits')
    g.save(record/'software.json',g.core.software())
    g.save(record/'retry-policy.json',dict(seeds=list(SEEDS),max_iterations=MAXIT,epsilon=0.00005,nonlinearity='pow3',approach='concat',
        change='Explicit external random starts and higher iteration budget only; first converged valid fit in fixed seed order.',
        fallback_nonlinearity=False,selection_by_task_effect=False,completed_fits_untouched=True,
        concurrent_queue_slots=jobs,memory_reservations_GiB=dict(fixed=2,automatic=16,headroom=8)))
    selected={}
    # Prefer a previously selected retry; otherwise retain the original canonical fit.
    chosen_file=root/'selected_ica_fits.json'
    earlier=read(chosen_file).get('fits',{}) if chosen_file.is_file() else {}
    pending=[]
    for plan in g.plans():
        mode,roi,dim=plan;key=f'{mode}/{roi}/dim-{dim}'
        marker=root/'ica'/mode/roi/f'dim-{dim}/complete.json'
        result=earlier.get(key) or (read(marker) if marker.is_file() else None)
        if result:
            selected[key]=validate_completed(result,root,verified);state['fits'][key]=selected[key]
        else:
            state['fits'][key]=dict(state='queued',mode=mode,roi=roi,dimensionality=dim)
            pending.append(plan)
    g.save(record/'state.json',state)
    print(f'Reusing {len(selected)} completed fits and all 490 prepared runs/1584 masked inputs; {len(pending)} fits queued.',flush=True)
    gate=g.inputs.MemoryGate();lock=threading.Lock()
    def save_selection():g.save(chosen_file,dict(preparation_source_record=str(prep_source),fits=selected))
    save_selection()
    state['state']='resuming_remaining_ICA_fits';g.save(record/'state.json',state)
    def fit(plan):
        mode,roi,dim=plan;key=f'{mode}/{roi}/dim-{dim}'
        def on_attempt(attempt):
            with lock:
                state['fits'][key]=dict(state='running',**{k:attempt[k] for k in ('mode','roi','dimensionality','seed')})
                g.save(record/'state.json',state)
        try:
            with gate.reserve(16 if dim==0 else 2):
                result=retry_fit(plan,root,masks,slices,reference,record,on_attempt)
        except Exception as exc:
            result=dict(state='failed',mode=mode,roi=roi,dimensionality=dim,error=str(exc))
            (record/'errors'/(key.replace('/','_')+'.txt')).write_text(traceback.format_exc())
        with lock:
            state['fits'][key]=result
            if result['state']=='complete':selected[key]=result
            save_selection();g.save(record/'state.json',state)
        print('ICA '+key+': '+result['state'],flush=True)
    # Automatic fits run independently even if a fixed-dimensionality retry fails.
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for future in as_completed([pool.submit(fit,p) for p in pending]):future.result()
    # A consistent, paired normalization pilot is kept separate from the 24 baseline fits.
    state['state']='running_native_variance_normalization_pilot';g.save(record/'state.json',state)
    from group_ica_vn_pilot import run_pilot
    run_pilot(root,record,state,masks,slices,reference,selected)
    g.core.require(len(selected)==24,'Some ICA fits remain incomplete; all completed fits and retry logs preserved')
    state['state']='generating_component_review';g.save(record/'state.json',state)
    g.review_outputs(list(selected.values()),record,masks,reference,state)
    if any(v['state']!='complete' for v in state['native_variance_normalization_pilot'].values()):
        state['state']='paired_ICA_complete_VN_pilot_incomplete';g.save(record/'state.json',state)
    (record/'README.md').write_text((record/'README.md').read_text()+'\n\nNative MELODIC variance-normalization right-insula dim-5 pilot: see `native-variance-normalization-pilot.json`. '
        'This is a separate paired diagnostic, not a change to selected baseline fits.\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    parser.add_argument('--jobs',type=int,default=11)
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args();g.core.require(1<=args.jobs<=11,'Resume supports 1–11 concurrent queue slots')
    records=g.core.ROOT/'code/insula_ica/records'
    candidates=sorted(p for p in records.glob('linux2-paired-group-ica-*') if (p/'state.json').is_file())
    g.core.require(args.record is not None or candidates,'No prepared paired ICA record found')
    source=(args.record or candidates[-1]).resolve()
    previous=read(source/'state.json');root=Path(previous['output_root']).resolve()
    g.core.require(root.is_relative_to((g.core.ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('paired-group-ica-'),'Invalid prepared analysis root')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    lockpath=g.core.ROOT/'code/insula_ica/output/paired-group-ica.lock';lockpath.parent.mkdir(parents=True,exist_ok=True)
    with lockpath.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        record=records/('linux2-paired-group-ica-'+stamp+'-resume');record.mkdir()
        for name in ('errors','figures','attempts'):(record/name).mkdir()
        state=dict(state='checking_prepared_inputs',runs=previous.get('runs',{}),fits={},output_root=str(root),source_record=str(source),
            started_utc=datetime.now(timezone.utc).isoformat(),code_commit=g.subprocess.check_output(['git','rev-parse','HEAD'],cwd=g.core.ROOT,text=True).strip())
        g.save(record/'state.json',state)
        try:run(source,root,record,state,args.jobs)
        except Exception as exc:
            state.update(state='failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        g.save(record/'state.json',state)
        if args.publish:
            try:g.publish(record,'Resume unfinished insula ICA with explicit deterministic optimization retries')
            except Exception as exc:
                (record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
        print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state']=='paired_ICA_complete_review_pending' else 1)


if __name__=='__main__':main()
