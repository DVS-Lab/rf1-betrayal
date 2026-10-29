#!/usr/bin/env python3
"""Corrected sm-0 left/right/bilateral paired-scale ICA from 490 completed FEATs."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import traceback
import uuid

import nibabel as nib
import numpy as np

import run_group_ica as g
import resume_group_ica as retry
from run_bold_controls import save, publish


ROOT=g.core.ROOT
RECORDS=ROOT/'code/insula_ica/records'


def latest_source():
    candidates=sorted(p for p in RECORDS.glob('linux2-ugr-TR-repair-*') if (p/'state.json').is_file())
    g.core.require(bool(candidates),'No corrected timing-repair record found')
    repair=candidates[-1];state=json.loads((repair/'state.json').read_text())
    g.core.require(state['state']=='TR_repair_complete_QC_review_pending','Latest timing correction did not complete')
    return Path(state['combined_sm0_record']),repair


def select_inputs(source):
    state=json.loads((source/'state.json').read_text())
    g.core.require(state['state']=='sm0_FEAT_complete_QC_pending' and len(state['runs'])==490 and
                   all(v['state']=='complete' for v in state['runs'].values()),'All 490 corrected sm-0 FEATs must be complete')
    base=Path(state['source_record']);original=json.loads((base/'accepted-inputs.json').read_text())
    root=Path(state['output_root']);rows=[]
    for row in original:
        feat=root/'feat'/(row['key']+'_sm-0.feat');marker=json.loads((feat/'sm0_complete.json').read_text())
        qc=json.loads((source/'qc'/(row['key']+'.json')).read_text())
        g.core.require(marker['qc']==qc and qc['finite_residuals'] and qc['exact_original_nuisance_matrix'],'Incomplete/unverified sm-0 QC')
        g.core.require(all(g.core.identity(p['path'],p['sha256'] is not None)==p for p in marker['outputs']),'Completed sm-0 output changed: '+row['key'])
        settings=g.core.fsl_settings(feat/'design.fsf');tr=float(settings['fmri(tr)'])
        g.core.require(float(settings['fmri(smooth)'])==0,'Source smoothing is not zero')
        corrected=row['key']=='sub-10817_ses-01_task-ugr_run-1'
        g.core.require(np.isclose(tr,1.710 if corrected else 1.615,rtol=0,atol=1e-5),'Unexpected run-specific TR')
        if corrected:g.core.require(qc.get('task_design_regenerated_for_TR'),'TR-corrected task design not verified')
        else:g.core.require(qc.get('complete_design_preserved'),'Other source model unexpectedly changed')
        rows.append(dict(row,manuscript_original_feat=row['original_feat'],original_feat=str(feat),
            active_feat=str(feat),source_fsf=str(feat/'design.fsf'),filtered_data=str(feat/'filtered_func_data.nii.gz'),
            native_residual=str(feat/'stats/res4d.nii.gz'),tr=tr,smoothing_mm=0,
            source_qc=str(source/'qc'/(row['key']+'.json')),source_selection_record=str(source),
            model_reference='Validated completed sm-0 full model; corrected task sampling for authorized 1.710-second UGR case.'))
    cohort=set((ROOT/'code/sublist_n132.txt').read_text().split())
    g.core.require(len(rows)==len({r['key'] for r in rows})==490 and len(cohort)==132,'Unexpected cohort/run identities')
    g.core.require(Counter(r['task'] for r in rows)==dict(trust=241,ugr=249),'Task run counts changed')
    for task in ('trust','ugr'):g.core.require({r['subject'] for r in rows if r['task']==task}==cohort,'Task cohort changed')
    return rows


def render_partial(results,record,masks,reference):
    for result in results:
        if result['state']!='complete':continue
        maps=nib.load(result['oriented_maps']).get_fdata();dim=result['dimensionality'];roi=result['roi']
        for start in range(0,maps.shape[3],12):
            path=record/'figures'/f'{result["mode"]}_{roi}_dim-{dim}_page-{start//12+1}.png'
            g.component_figure(maps[...,start:start+12],masks[roi],reference.get_fdata(),reference.affine,path,
                f'Unsmoothed, native VN | {result["mode"]} | {roi} | requested {dim}',first_component=start+1)


def run(args,root,record,state):
    rows=select_inputs(args.source_record);reference=nib.load(Path(rows[0]['active_feat'])/'mean_func.nii.gz')
    save(record/'accepted-inputs.json',rows);save(record/'software.json',g.core.software())
    policy=dict(cohort_n=132,runs=490,trust_runs=241,ugr_runs=249,smoothing_mm=0,
        residual_modes=list(g.inputs.MODES),rois=list(g.core.ROIS),dimensions=[3,5,7,0],native_MELODIC_variance_normalization=True,
        global_input_amplitude_multiplier=1000,per_run_demean_SD_and_equal_run_energy=True,
        no_temporal_resampling=True,acquisition_TRs_preserved=True,common_TR_override=False,
        mixed_TR_frequency_report_interpretable=False,native_MELODIC_web_report=False,
        temporal_diagnostics='Sample-lag ACF only; no common-Hz interpretation.',
        source_record=str(args.source_record),corrected_run='sub-10817_ses-01_task-ugr_run-1',
        new_regression=False,new_temporal_filter=False,new_FEAT=False,new_FILM=False,new_fmriprep=False,
        seeds=list(retry.SEEDS),iterations_per_start=retry.MAXIT,nonlinearity='pow3',epsilon=.00005,
        prepare_slots=args.jobs,fixed_ICA_slots=18,auto_ICA_slots=1,auto_memory_reservation_GiB=32,
        fixed_threads=4,automatic_threads=8,automatic_PCA='Classic full group PCA, serial automatic fits; no MIGP substitution.')
    frozen=root/'policy.json'
    if frozen.is_file():
        old=json.loads(frozen.read_text());old['prepare_slots']=policy['prepare_slots']
        g.core.require(old==policy,'Unsmoothed ICA policy changed; use a fresh output root')
    else:save(frozen,policy)
    save(record/'policy.json',policy)
    fsldir=Path(os.environ['FSLDIR']);atlas=fsldir/'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'
    xml=fsldir/'data/atlases/HarvardOxford-Cortical.xml';anatomical,counts,label=g.core.anatomical_masks(atlas,xml,reference)
    save(record/'atlas.json',dict(label=label,counts=counts,atlas=g.core.identity(atlas),xml=g.core.identity(xml)))
    state['state']='validating_corrected_sm0_sources';save(record/'state.json',state)
    for row in rows:
        g.core.require(Path(row['native_residual']).is_file(),'Completed native residual missing; no FILM will be launched')
        g.inputs.inspect_source(row,reference)
    g.core.require(shutil.disk_usage(root).free>=16*1024**3,'Need 16 GiB free for retained ICA preparation/outputs')
    # Verify the actual native-FILM decomposition on both tasks and the corrected
    # 1.710-second case before any group ICA is launched.
    state['state']='paired_scale_pilot';save(record/'state.json',state)
    pilot=[rows[0],next(r for r in rows if r['task']=='ugr'),next(r for r in rows if np.isclose(r['tr'],1.710))]
    checks=[]
    for row in pilot:
        feat=Path(row['active_feat']);qc=g.residuals.residual_qc(feat,feat,dict(row,prewhiten='1'),
            dict(output_root=str(root),residual_projection_tolerance=1e-4))
        g.core.require(qc['scale']=='voxelwise_prewhitened' and qc.get('max_whitened_decomposition_relative_error',1)<1e-4,'Native residual scale pilot failed')
        checks.append(qc)
    save(record/'scale-pilot.json',dict(passed=True,runs=checks))
    gate=g.inputs.MemoryGate();lock=threading.Lock()
    state['state']='preparing_unsmoothed_paired_residuals';save(record/'state.json',state)
    def prepare(row):
        try:
            with gate.reserve(1):
                g.core.require(Path(row['native_residual']).is_file(),'Native source disappeared; no FILM rerun permitted')
                qc=g.inputs.prepare_run(row,root,record,anatomical['bilateral'],reference)
            g.core.require(not qc['native_recovered_by_film'],'Unexpected FILM recovery')
            result=dict(state='complete')
        except Exception as exc:
            result=dict(state='failed',error=str(exc));(record/'errors'/(row['key']+'.txt')).write_text(traceback.format_exc())
        with lock:
            state['runs'][row['key']]=result;save(record/'state.json',state)
            if len(state['runs'])%25==0 or result['state']=='failed':print(f'Unsmoothed residuals {len(state["runs"])}/490 | {row["key"]}: {result["state"]}',flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for done in as_completed([pool.submit(prepare,row) for row in rows]):done.result()
    g.core.require(all(r['state']=='complete' for r in state['runs'].values()),'Residual preparation failed; no ICA launched')
    masks=g.fixed_masks(rows,root,record,anatomical,reference)
    state['state']='building_physically_masked_TR_preserving_inputs';save(record/'state.json',state)
    slices,_=g.sessions(rows,root,record,masks,anatomical['bilateral'],reference,amplitude=1000,preserve_timing=True)
    save(record/'timing-policy.json',dict(acquisition_tr_counts=dict(Counter(str(r['tr']) for r in rows)),
        frequency_spectra_interpretable=False,interpolation=False,physical_headers_preserved=True,
        note='Spatial temporal-concatenation ICA only. No common-TR time/frequency inference; no MELODIC web-report plots are generated.'))
    results=[]
    def fit(plan):
        mode,roi,dim=plan;key=f'{mode}/{roi}/dim-{dim}'
        try:
            with gate.reserve(32 if dim==0 else 2):
                result=retry.retry_fit(plan,root,masks,slices,reference,record,
                    variance_normalization=True,preserve_timing=True,threads=8 if dim==0 else 4)
        except Exception as exc:
            result=dict(mode=mode,roi=roi,dimensionality=dim,state='failed',error=str(exc))
            (record/'errors'/(key.replace('/','_')+'.txt')).write_text(traceback.format_exc())
        with lock:
            results.append(result);state['fits'][key]=result;save(record/'state.json',state)
        print('Unsmoothed ICA '+key+': '+result['state'],flush=True)
    state['state']='running_18_unsmoothed_fixed_ICA_fits';save(record/'state.json',state)
    with ThreadPoolExecutor(max_workers=18) as pool:
        for done in as_completed([pool.submit(fit,p) for p in g.plans() if p[2]!=0]):done.result()
    # The automatic fits run even if a fixed fit failed. They are serialized to
    # avoid the previous six simultaneous classic-PCA memory peaks.
    state['state']='running_6_serial_automatic_ICA_fits';save(record/'state.json',state)
    if Path('/proc/meminfo').is_file():
        total=next(int(x.split()[1])/1024**2 for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemTotal:'))
        g.core.require(total>=40,'Automatic PCA needs 32 GiB plus 8 GiB headroom; fixed fits preserved')
    for plan in g.plans():
        if plan[2]==0:fit(plan)
    state['state']='rendering_unsmoothed_ICA_review';save(record/'state.json',state)
    if all(r['state']=='complete' for r in results):
        g.review_outputs(results,record,masks,reference,state)
        state['state']='unsmoothed_ICA_complete_review_pending'
    else:
        render_partial(results,record,masks,reference);state['state']='unsmoothed_ICA_incomplete_review_errors'
    save(record/'state.json',state)
    with (record/'README.md').open('a') as stream:
        stream.write('\n\n## Corrected unsmoothed stream\n\n'
            'All source FEATs use zero smoothing. This includes the regenerated 1.710-second UGR task design. '
            'Both residual scales use native MELODIC variance normalization, in addition to the documented runwise balancing; '
            'inputs have a uniform ×1000 multiplier to avoid the absolute native-VN SD floor. '
            'TRs are preserved in each participant/task block. There is no temporal interpolation, common-TR override, new temporal filter or nuisance regression. '
            'MELODIC web reports are suppressed; review images and comparisons concern spatial maps. '
            'Mixed-TR concatenated frequency plots are not interpreted. Automatic fits use classic PCA, one at a time, 8 threads and a 32 GiB RAM reservation. '
            'Existing smoothed/baseline fits and first-level outputs were not overwritten. '
            f'State: {state["state"]}. See state.json for all 24 outcomes; failure logs/partial successful figures are retained.\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-record',type=Path);parser.add_argument('--resume-record',type=Path)
    parser.add_argument('--jobs',type=int,default=64);parser.add_argument('--publish',action='store_true');args=parser.parse_args()
    g.core.require(1<=args.jobs<=64,'Expected 1–64 preparation workers')
    if args.resume_record is None and args.source_record is None:
        previous=sorted(p for p in RECORDS.glob('linux2-unsmoothed-group-ica-*') if (p/'state.json').is_file())
        if previous:args.resume_record=previous[-1]
    if args.resume_record:
        old=json.loads((args.resume_record/'state.json').read_text());root=Path(old['output_root']);args.source_record=Path(old['source_record'])
    else:
        root=None
        if args.source_record is None:args.source_record,_=latest_source()
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    if root is None:root=ROOT/'derivatives/insula_ica'/('unsmoothed-group-ica-'+stamp)
    g.core.require(root.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('unsmoothed-group-ica-'),'Invalid unsmoothed ICA root')
    output=ROOT/'code/insula_ica/output';output.mkdir(exist_ok=True)
    with (output/'unsmoothed-group-ica.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        record=RECORDS/('linux2-unsmoothed-group-ica-'+stamp);record.mkdir()
        for name in ('runs','figures','attempts','errors'):(record/name).mkdir()
        root.mkdir(parents=True,exist_ok=True)
        state=dict(state='preflight',source_record=str(args.source_record),output_root=str(root),runs={},fits={},
            code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
        save(record/'state.json',state)
        try:run(args,root,record,state)
        except Exception as exc:
            state.update(state='failed',error=str(exc));save(record/'state.json',state)
            (record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        if args.publish:
            try:publish(record,'Run corrected unsmoothed paired-scale left/right/bilateral insula ICA')
            except Exception as exc:
                (record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
        print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state']=='unsmoothed_ICA_complete_review_pending' else 1)


if __name__=='__main__':main()
