#!/usr/bin/env python3
"""Paired BOLD-scale/native-FILM insula ICA, with accepted Trust replacements."""
import argparse
from collections import defaultdict
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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy.linalg import subspace_angles
from scipy.optimize import linear_sum_assignment

import insula_ica as core
import feat_residuals as residuals
import group_ica_inputs as inputs
from run_bold_controls import save, publish


DIMENSIONS=(5,3,7,0)
SEED=20261028
AUDIT=core.ROOT/'code/insula_ica/output/linux2-audit-20261026T213907Z-dc0e41/audit'


def plans():
    return [(mode,roi,dim) for mode in inputs.MODES for roi in core.ROIS for dim in DIMENSIONS]


def compare_maps(first,second,*,labels=('joint_bold','native_film')):
    """Match signed ICA maps by absolute spatial correlation, retaining all pairs."""
    first=np.asarray(first,dtype=float);second=np.asarray(second,dtype=float)
    core.require(first.ndim==second.ndim==2 and first.shape[0]==second.shape[0],'Map comparison grid mismatch')
    a=first-first.mean(axis=0);b=second-second.mean(axis=0)
    core.require((np.linalg.norm(a,axis=0)>0).all() and (np.linalg.norm(b,axis=0)>0).all(),'Constant component map')
    corr=(a.T@b)/np.outer(np.linalg.norm(a,axis=0),np.linalg.norm(b,axis=0))
    left,right=linear_sum_assignment(-np.abs(corr))
    result=dict(joint_bold_components=first.shape[1],native_film_components=second.shape[1],
        unmatched_joint_bold_components=[k+1 for k in range(first.shape[1]) if k not in left],
        unmatched_native_film_components=[k+1 for k in range(second.shape[1]) if k not in right],
        all_pair_spatial_r=corr.tolist(),matches=[dict(joint_bold_component=int(i+1),
        native_film_component=int(j+1),spatial_r=float(corr[i,j]),absolute_spatial_r=float(abs(corr[i,j])),
        native_sign_for_joint_comparison=1 if corr[i,j]>=0 else -1) for i,j in zip(left,right)],
        centered_subspace_angles_degrees=np.degrees(subspace_angles(a,b)).tolist(),
        matching_rule='One-to-one maximum total absolute spatial correlation; all pairwise correlations also retained.')
    if labels!=('joint_bold','native_film'):
        for old,new in zip(('joint_bold','native_film'),labels):
            for prefix in ('','unmatched_'):
                result[prefix+new+'_components']=result.pop(prefix+old+'_components')
        for match in result['matches']:
            match[labels[0]+'_component']=match.pop('joint_bold_component')
            match[labels[1]+'_component']=match.pop('native_film_component')
            match[labels[1]+'_sign_for_'+labels[0]+'_comparison']=match.pop('native_sign_for_joint_comparison')
    return result


def write_nifti(path,data,affine,tr=1.615):
    image=nib.Nifti1Image(np.asarray(data),affine)
    image.set_qform(affine,1);image.set_sform(affine,1)
    image.header.set_xyzt_units('mm','sec')
    if data.ndim==4:image.header.set_zooms((*image.header.get_zooms()[:3],tr))
    nib.save(image,path)


def freeze_nifti(path,data,affine):
    if path.exists():
        old=nib.load(path)
        core.require(old.shape==data.shape and np.allclose(old.affine,affine,atol=1e-5) and
                     np.array_equal(np.asanyarray(old.dataobj),data),'Previously frozen image changed: '+str(path))
    else:write_nifti(path,data,affine)


def validate_masked_4d(path,mask_image,frames):
    """Read the written input: every voxel outside this specific ROI must be zero."""
    image=nib.load(path);mask=mask_image.get_fdata()>0
    core.require(core.same_grid(image,mask_image) and image.shape==(*mask.shape,frames),'Masked 4D input grid/timing mismatch')
    core.require(np.isclose(image.header.get_zooms()[3],1.615,atol=1e-5),'Masked 4D input TR mismatch')
    data=np.asanyarray(image.dataobj)
    core.require(np.isfinite(data).all(),'Nonfinite masked 4D input')
    core.require(np.count_nonzero(data[~mask])==0,'Nonzero data outside the specific insula ROI: '+str(path))
    core.require((data[mask].std(axis=1)>1e-8).all(),'Constant data inside the specific insula ROI')
    return dict(physical_4d_mask_verified=True,outside_roi_nonzero_values=0,frames=frames,
                inside_roi_voxels=int(mask.sum()),output=core.identity(path))


def verify_ica_inputs(mode,roi,root):
    folder=root/'sessions'/mode/roi
    manifest=json.loads((folder/'verified-inputs.json').read_text())
    listing=(folder/'inputs.txt').read_text().splitlines()
    core.require(manifest['mode']==mode and manifest['roi']==roi and len(listing)==264,'Wrong ICA input collection')
    core.require(manifest['mask_sha256']==core.sha(root/'masks'/(roi+'_crop.nii.gz')),'Physical input mask changed')
    core.require(listing==[v['output']['path'] for v in manifest['inputs']],'ICA input listing changed')
    for item in manifest['inputs']:
        core.require(item['physical_4d_mask_verified'] and item['outside_roi_nonzero_values']==0,'Unverified physical ROI masking')
        core.require(core.identity(item['output']['path'])==item['output'],'Previously verified 4D input changed')
    return core.sha(folder/'verified-inputs.json')


def pilot_checks(rows,root,record,reference):
    trust=[r for r in rows if r['task']=='trust'];ugr=[r for r in rows if r['task']=='ugr']
    chosen=[trust[0],next(r for r in trust if r['key']=='sub-10834_ses-01_task-trust_run-1'),
            next(r for r in trust if r['key']=='sub-10720_ses-01_task-trust_run-2'),ugr[0],ugr[-1]]
    results=[]
    for row in chosen:
        inputs.inspect_source(row,reference)
        feat,replayed=inputs.ensure_native(row,root,require_autocorr=True)
        qc=residuals.residual_qc(feat,feat,dict(row,prewhiten='1'),
            dict(output_root=str(root),residual_projection_tolerance=1e-4))
        core.require(qc['scale']=='voxelwise_prewhitened' and
            qc.get('max_whitened_decomposition_relative_error',1)<1e-4,'Pilot native FILM decomposition failed')
        qc.update(native_feat=str(feat),native_replayed=replayed)
        results.append(qc)
        print('Paired-scale pilot passed:',row['key'],flush=True)
    save(record/'scale-pilot.json',dict(passed=True,runs=results,
        joint_bold_definition='Y minus its temporal mean minus the full saved design times the existing joint FILM PEs.',
        native_film_definition='Native FILM voxelwise-prewhitened innovations, verified using the matching voxelwise whitened design.',
        extra_regression_or_temporal_filter=False))


def fixed_masks(rows,root,record,anatomical,reference):
    union=anatomical['bilateral'];common=np.ones(int(union.sum()),dtype=bool)
    coverage=[]
    for row in rows:
        valid=np.load(root/'run_arrays'/(row['key']+'_valid.npy'))
        common &= valid
        coverage.append(dict(key=row['key'],**{roi:float(valid[anatomical[roi][union]].mean()) for roi in core.ROIS}))
    full_valid=np.zeros(union.shape,dtype=bool);full_valid[union]=common
    masks={roi:mask & full_valid for roi,mask in anatomical.items()}
    summary={roi:dict(anatomical_voxels=int(anatomical[roi].sum()),common_valid_voxels=int(mask.sum()),
        retained_fraction=float(mask.sum()/anatomical[roi].sum())) for roi,mask in masks.items()}
    save(record/'coverage.json',dict(per_run=coverage,common=summary,
        definition='Same fitted-mask/finite/nonconstant support across every accepted run and both residual representations.'))
    core.require(all(s['retained_fraction']>=.8 for s in summary.values()),'Common insula coverage below the existing 80% gate; inspect coverage.json, no automatic exclusions')
    destination=root/'masks';destination.mkdir(exist_ok=True)
    for roi in core.ROIS:
        for suffix,data in (('',masks[roi]),('_lost',anatomical[roi]&~masks[roi])):
            path=destination/(roi+suffix+'.nii.gz')
            if path.exists():
                previous=nib.load(path)
                core.require(core.same_grid(previous,reference) and np.array_equal(previous.get_fdata()>0,data),'Previously frozen mask changed')
            else:core.save_image(path,data.astype('uint8'),reference)
    return masks


def sessions(rows,root,record,masks,anatomical_union,reference,*,modes=inputs.MODES,amplitude=1.0):
    slices,affine=inputs.crop_geometry(reference,masks['bilateral'])
    shape=masks['bilateral'][slices].shape
    groups=defaultdict(list)
    for row in rows:groups[(row['subject'],row['task'])].append(row)
    core.require(len(groups)==264,'Expected 264 participant/task PCA input blocks')
    kept=masks['bilateral'][anatomical_union]
    for roi in core.ROIS:
        freeze_nifti(root/'masks'/(roi+'_crop.nii.gz'),masks[roi][slices].astype('uint8'),affine)
    ref_data=reference.get_fdata()[slices]
    freeze_nifti(root/'reference_crop.nii.gz',ref_data.astype('float32'),affine)
    provenance=[];verification=[]
    for mode in modes:
        folders={roi:root/'sessions'/mode/roi for roi in core.ROIS}
        listings={roi:[] for roi in core.ROIS};validated={roi:[] for roi in core.ROIS}
        for folder in folders.values():folder.mkdir(parents=True,exist_ok=True)
        for (subject,task),rr in sorted(groups.items()):
            rr=sorted(rr,key=lambda r:r['run'])
            arrays=[np.load(root/'run_arrays'/(r['key']+'_'+mode+'.npy'))[kept] for r in rr]
            values,weights=inputs.balanced_session(arrays)
            for roi,folder in folders.items():
                destination=folder/f'sub-{subject}_task-{task}.nii.gz'
                marker=destination.with_name(destination.name+'.complete.json')
                signature=dict(run_arrays=[core.sha(root/'run_arrays'/(r['key']+'_complete.json')) for r in rr],
                    mask_sha256=core.sha(root/'masks'/(roi+'_crop.nii.gz')),mode=mode,roi=roi)
                if amplitude!=1:signature['global_amplitude_multiplier']=amplitude
                if marker.exists():
                    old=json.loads(marker.read_text());check=old['verification']
                    core.require(old['signature']==signature and check['physical_4d_mask_verified'] and
                        core.identity(destination)==check['output'],'Balanced ICA session changed')
                else:
                    core.require(not destination.exists(),'Partial balanced input exists; preserve it')
                    data=np.zeros((*shape,values.shape[1]),dtype='float32')
                    data[masks[roi][slices]]=values[masks[roi][masks['bilateral']]]*amplitude
                    write_nifti(destination,data,affine)
                    check=validate_masked_4d(destination,nib.load(root/'masks'/(roi+'_crop.nii.gz')),values.shape[1])
                    save(marker,dict(signature=signature,verification=check))
                listings[roi].append(str(destination.resolve()));validated[roi].append(check)
            provenance.append(dict(mode=mode,subject=subject,task=task,runs=[r['key'] for r in rr],
                frames=[r['nvol'] for r in rr],within_session_weights=weights))
        for roi,folder in folders.items():
            (folder/'inputs.txt').write_text('\n'.join(listings[roi])+'\n')
            manifest=dict(mode=mode,roi=roi,mask_sha256=core.sha(root/'masks'/(roi+'_crop.nii.gz')),inputs=validated[roi])
            save(folder/'verified-inputs.json',manifest)
            verification.append(dict(mode=mode,roi=roi,physically_masked_4d_inputs=len(validated[roi]),
                outside_roi_nonzero_values=0,verified_inputs_sha256=verify_ica_inputs(mode,roi,root)))
    geometry=dict(full_shape=list(reference.shape),full_affine=reference.affine.tolist(),
        crop_start=[s.start for s in slices],crop_stop=[s.stop for s in slices],crop_affine=affine.tolist(),
        policy='Array crop only. MNI voxel coordinates unchanged; no resampling. Component maps are also saved on the full original grid.')
    save(record/'crop.json',geometry);save(record/'contributions.json',provenance)
    save(record/'physical-4d-masking.json',verification)
    return slices,affine


def component_figure(maps,mask,background,affine,path,label,first_component=1):
    n=maps.shape[3];fig,axes=plt.subplots(n,3,figsize=(9,2.5*n),squeeze=False,facecolor='black')
    inside=background[np.isfinite(background)&(background>0)]
    low,high=np.percentile(inside,[1,99])
    for k in range(n):
        values=maps[...,k];indices=np.argwhere(mask)
        peak=indices[np.argmax(np.abs(values[mask]))];vmax=float(np.max(np.abs(values[mask])))
        xyz=nib.affines.apply_affine(affine,peak)
        display=np.ma.masked_where(np.abs(values)<.3*vmax,values)
        for axis in range(3):
            sl=[slice(None)]*3;sl[axis]=int(peak[axis]);sl=tuple(sl)
            ax=axes[k,axis];ax.imshow(background[sl].T,origin='lower',cmap='gray',vmin=low,vmax=high)
            ax.imshow(display[sl].T,origin='lower',cmap='coolwarm',vmin=-vmax,vmax=vmax,alpha=.85)
            ax.set_title(f'IC {k+first_component}; {"xyz"[axis]}={xyz[axis]:.1f} mm',color='white',fontsize=9);ax.axis('off')
    fig.suptitle(label+'\nLargest absolute voxel defines sign. Display: |loading| ≥ 30% of peak, no inference threshold.',color='white',fontsize=10)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(path,dpi=100,facecolor='black');plt.close(fig)


def melodic_command(mode,roi,dim,root,*,output_folder=None,seed=SEED,maxit=2000,variance_normalization=False):
    argv=['melodic','-i',str(root/'sessions'/mode/roi/'inputs.txt'),'-o',str(output_folder or root/'ica'/mode/roi/f'dim-{dim}'),
        '--tr=1.615','--nobet','-a','concat','-m',str(root/'masks'/(roi+'_crop.nii.gz')),
        '--bgimage='+str(root/'reference_crop.nii.gz'),'--report','--Oall','-d',str(dim),
        '--disableMigp',*([] if variance_normalization else ['--vn']),'--update_mask','--seed='+str(seed),'--maxit='+str(maxit)]
    if output_folder is None and seed==SEED and maxit==2000 and not variance_normalization:argv.append('--maxrestart=10')
    else:argv.extend(['--eps=0.00005','--nl=pow3'])
    return argv


def fit_one(mode,roi,dim,root,masks,slices,reference,*,output_folder=None,seed=SEED,maxit=2000,variance_normalization=False):
    verified_inputs_sha256=verify_ica_inputs(mode,roi,root)
    argv=melodic_command(mode,roi,dim,root,output_folder=output_folder,seed=seed,maxit=maxit,variance_normalization=variance_normalization)
    folder=Path(argv[argv.index('-o')+1]);folder.parent.mkdir(parents=True,exist_ok=True)
    marker=folder/'complete.json'
    if marker.is_file():
        previous=json.loads(marker.read_text())
        core.require(previous['argv']==argv and previous['verified_inputs_sha256']==verified_inputs_sha256 and previous['mask_sha256']==core.sha(root/'masks'/(roi+'.nii.gz')) and
            previous['map_identity']==core.identity(previous['oriented_maps']),'Completed ICA inputs/maps changed')
        return dict(previous,reused_completed_ica=True)
    core.require(not folder.exists(),'Partial/completed ICA directory already exists; preserve it')
    env=dict(os.environ,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
    log=folder.parent/(folder.name+'.log')
    with log.open('w') as stream:subprocess.run(argv,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
    core.require((folder/'melodic_IC.nii.gz').is_file(),'MELODIC produced no component maps (possible nonconvergence); inspect '+str(log))
    mask_image=nib.load(root/'masks'/(roi+'_crop.nii.gz'));mask=mask_image.get_fdata()>0
    actual=nib.load(folder/'mask.nii.gz')
    core.require(core.same_grid(actual,mask_image) and np.array_equal(actual.get_fdata()>0,mask),'MELODIC changed fixed mask')
    image=nib.load(folder/'melodic_IC.nii.gz');data=image.get_fdata()
    if data.ndim==3:data=data[...,None]
    core.require(core.same_grid(image,mask_image) and data.ndim==4 and data.shape[:3]==mask.shape and
        (dim==0 or data.shape[3]==dim),'Unexpected ICA dimensions/grid')
    actual_dim=data.shape[3]
    core.require(np.isfinite(data).all(),'Nonfinite ICA maps')
    spatial=data[mask]-data[mask].mean(axis=0)
    core.require(np.linalg.matrix_rank(np.column_stack([np.ones(mask.sum()),spatial]))==actual_dim+1,'Nonestimable spatial ICA basis')
    data,signs=core.orient_maps(data,mask)
    write_nifti(folder/'oriented_IC_crop.nii.gz',data.astype('float32'),image.affine)
    full=np.zeros((*reference.shape,actual_dim),dtype='float32');full[slices]=data
    core.save_image(folder/'oriented_IC.nii.gz',full,reference)
    result=dict(mode=mode,roi=roi,dimensionality=dim,actual_components=actual_dim,
        dimensionality_selection='MELODIC automatic estimation' if dim==0 else 'prespecified',argv=argv,seed=seed,
        max_iterations=maxit,convergence_epsilon=0.00005,nonlinearity='pow3',variance_normalization=variance_normalization,threads_per_fit=4,
        oriented_maps=str(folder/'oriented_IC.nii.gz'),mask_sha256=core.sha(root/'masks'/(roi+'.nii.gz')),
        signs=signs,map_identity=core.identity(folder/'oriented_IC.nii.gz'),verified_inputs_sha256=verified_inputs_sha256,
        no_new_temporal_filter=True,no_new_nuisance_regression=True)
    save(folder/'complete.json',result)
    return result


def run(args,root,record,state):
    software=core.software();save(record/'software.json',software)
    help_result=subprocess.run(['melodic','--help'],capture_output=True,text=True)
    help_text=help_result.stdout+help_result.stderr;(record/'melodic-help.txt').write_text(help_text)
    for option in ('--update_mask','--disableMigp','--vn'):
        core.require(option in help_text,'Required installed MELODIC option unavailable: '+option)
    rows=inputs.select_inputs(args.audit)
    save(record/'accepted-inputs.json',rows)
    policy=dict(analysis_policy_version=1,cohort_n=132,trust_runs=241,ugr_runs=249,
        residual_modes=list(inputs.MODES),rois=list(core.ROIS),primary_dimensionality=5,sensitivity_dimensionalities=[3,7],
        automatic_dimensionality=True,automatic_dimension_code=0,auto_fit_ram_reservation_GiB=16,
        user_accepted_flagged_runs=sorted(inputs.ACCEPTED_FLAGS),both_representations_requested=True,
        physical_4d_masking='Separate saved left/right/bilateral 4D inputs; outside-ROI values verified exactly zero before MELODIC.',
        nuisance_regression='Existing complete joint task+nuisance fit only; no new regression.',
        interpretation=dict(joint_bold='Original filtered-BOLD scale joint GLS residual.',native_film='Voxelwise-prewhitened innovation scale; covariance differs.'),
        analysis_scope='Paired spatial ICA and temporal diagnostics only. Whole-brain connectivity/task inference remain separate.',
        native_residual_scale_review='User explicitly requested paired native-FILM and joint-BOLD ICA comparison; innovation-scale covariance is labeled separately.',
        prepare_ram_reservations_GiB=dict(existing_native=1,missing_native_FILM_replay=4,available_headroom=8),
        prepare_worker_limit=args.jobs,ica_parallel_limit=args.ica_jobs,automatic_ica_parallel_limit=args.auto_ica_jobs)
    frozen_policy=root/'policy.json'
    if frozen_policy.exists():
        old=json.loads(frozen_policy.read_text())
        for key in ('prepare_worker_limit','ica_parallel_limit','automatic_ica_parallel_limit'):old[key]=policy[key]
        core.require(old==policy,'Analysis policy changed; use a fresh analysis root')
    else:save(frozen_policy,policy)
    save(record/'policy.json',policy)
    reference=nib.load(Path(rows[0]['active_feat'])/'mean_func.nii.gz')
    fsldir=Path(os.environ['FSLDIR'])
    atlas=fsldir/'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'
    xml=fsldir/'data/atlases/HarvardOxford-Cortical.xml'
    anatomical,counts,label=core.anatomical_masks(atlas,xml,reference)
    save(record/'atlas.json',dict(label=label,counts=counts,atlas=core.identity(atlas),xml=core.identity(xml)))
    # Validate all retained sources before starting any new FILM or ICA job.
    state.update(state='checking_490_retained_fits',accepted_runs=len(rows));save(record/'state.json',state)
    inputs.verify_frozen_ugr(rows,args.audit)
    for index,row in enumerate(rows,1):
        try:inputs.inspect_source(row,reference)
        except Exception as exc:raise ValueError(row['key']+': '+str(exc)) from exc
        if index%50==0 or index==len(rows):print(f'Validated source models/confounds: {index}/490',flush=True)
    missing=[r for r in rows if not Path(r['native_residual']).is_file()]
    estimated_bytes=3*sum(int(np.prod(reference.shape))*r['nvol']*4 for r in missing)+8*1024**3
    available_bytes=shutil.disk_usage(root).free
    save(record/'disk-preflight.json',dict(available_bytes=available_bytes,
        conservative_required_free_bytes=estimated_bytes,missing_native_residuals=len(missing)))
    core.require(available_bytes>=estimated_bytes,'Insufficient free space for retained native FILM recovery outputs')
    state['state']='paired_scale_pilot';save(record/'state.json',state)
    pilot_checks(rows,root,record,reference)
    state['state']='preparing_paired_insula_residuals';save(record/'state.json',state)
    gate=inputs.MemoryGate();state_lock=threading.Lock()
    def prepare(row):
        try:
            with gate.reserve(1 if Path(row['native_residual']).is_file() else 4):
                qc=inputs.prepare_run(row,root,record,anatomical['bilateral'],reference)
            result=dict(state='complete',native_replayed=qc['native_recovered_by_film'])
        except Exception as exc:
            result=dict(state='failed',error=str(exc));(record/'errors'/(row['key']+'.txt')).write_text(traceback.format_exc())
        with state_lock:
            state['runs'][row['key']]=result;save(record/'state.json',state)
            done=len(state['runs'])
            if done%10==0 or result['state']=='failed':print(f'Paired inputs {done}/490 | {row["key"]}: {result["state"]}',flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(prepare,row) for row in rows]):future.result()
    failures=[key for key,r in state['runs'].items() if r['state']=='failed']
    core.require(not failures,'Input preparation failed for '+str(len(failures))+' runs; inspect errors, no ICA started')
    masks=fixed_masks(rows,root,record,anatomical,reference)
    state['state']='building_balanced_PCA_sessions';save(record/'state.json',state)
    slices,_=sessions(rows,root,record,masks,anatomical['bilateral'],reference)
    save(record/'ica-plan.json',[dict(mode=m,roi=r,dimensionality=d,argv=melodic_command(m,r,d,root)) for m,r,d in plans()])
    state['state']='running_18_fixed_dimensionality_ICA_fits';save(record/'state.json',state)
    results=[]
    def fit(plan):
        m,r,d=plan;key=f'{m}/{r}/dim-{d}'
        try:
            result=fit_one(m,r,d,root,masks,slices,reference);result['state']='complete'
        except Exception as exc:
            result=dict(mode=m,roi=r,dimensionality=d,state='failed',error=str(exc))
            (record/'errors'/(key.replace('/','_')+'.txt')).write_text(traceback.format_exc())
        with state_lock:
            results.append(result);state['fits'][key]=result;save(record/'state.json',state)
        print('ICA '+key+': '+result['state'],flush=True)
    with ThreadPoolExecutor(max_workers=args.ica_jobs) as pool:
        for future in as_completed([pool.submit(fit,p) for p in plans() if p[2]!=0]):future.result()
    core.require(all(r['state']=='complete' for r in results),'Some fixed-dimensionality ICA fits failed; completed fits remain preserved')
    state['state']='running_6_automatic_dimensionality_ICA_fits';save(record/'state.json',state)
    auto_gate=inputs.MemoryGate()
    def automatic_fit(plan):
        with auto_gate.reserve(16):fit(plan)
    with ThreadPoolExecutor(max_workers=args.auto_ica_jobs) as pool:
        for future in as_completed([pool.submit(automatic_fit,p) for p in plans() if p[2]==0]):future.result()
    core.require(all(r['state']=='complete' for r in results),'Some ICA fits failed; completed fits remain preserved')
    review_outputs(results,record,masks,reference,state)


def review_outputs(results,record,masks,reference,state):
    """Render and compare the selected completed fits, including preserved retry paths."""
    bykey={(r['mode'],r['roi'],r['dimensionality']):r for r in results}
    core.require(set(bykey)==set(plans()) and all(r['state']=='complete' for r in results),'All 24 selected fits must pass before final comparisons')
    # Figure rendering is serial: Matplotlib has global state.
    for result in results:
        maps=nib.load(result['oriented_maps']).get_fdata();roi=result['roi'];dim=result['dimensionality'];mode=result['mode']
        label=f'automatic: {maps.shape[3]} components' if dim==0 else f'{dim} components'
        result['figure_pages']=[]
        for start in range(0,maps.shape[3],12):
            suffix=f'_page-{start//12+1}' if maps.shape[3]>12 else ''
            path=record/'figures'/f'{mode}_{roi}_dim-{dim}{suffix}.png'
            component_figure(maps[...,start:start+12],masks[roi],reference.get_fdata(),reference.affine,
                path,f'{mode} | {roi} insula | {label}',first_component=start+1)
            result['figure_pages'].append(str(path.relative_to(record)))
    comparisons=[]
    for roi in core.ROIS:
        for dim in DIMENSIONS:
            a=nib.load(bykey['joint_bold',roi,dim]['oriented_maps']).get_fdata()[masks[roi]]
            b=nib.load(bykey['native_film',roi,dim]['oriented_maps']).get_fdata()[masks[roi]]
            comparisons.append(dict(roi=roi,dimensionality=dim,**compare_maps(a,b)))
    save(record/'between-scale-spatial-comparisons.json',comparisons)
    stability=[]
    for mode in inputs.MODES:
        for roi in core.ROIS:
            primary=nib.load(bykey[mode,roi,5]['oriented_maps']).get_fdata()[masks[roi]]
            for dim in (3,7,0):
                other=nib.load(bykey[mode,roi,dim]['oriented_maps']).get_fdata()[masks[roi]]
                corr=np.corrcoef(primary.T,other.T)[:5,5:]
                stability.append(dict(mode=mode,roi=roi,sensitivity_dim=dim,actual_components=other.shape[1],all_pair_spatial_r=corr.tolist()))
    save(record/'dimensionality-stability.json',stability)
    state.update(state='paired_ICA_complete_review_pending',finished_utc=datetime.now(timezone.utc).isoformat())
    save(record/'state.json',state)
    links='\n'.join(f'- [{r["mode"]}, {r["roi"]}, requested {r["dimensionality"]}, obtained {r["actual_components"]} components, page {j+1}]({path})'
        for r in sorted(results,key=lambda r:(r['mode'],r['roi'],r['dimensionality'])) for j,path in enumerate(r['figure_pages']))
    (record/'README.md').write_text('# Paired residual-scale insula ICA\n\nN=132; Trust 241 and UGR 249 runs; no new exclusions. '
        'Left/right/bilateral masks share common valid support across both residual representations and all runs. '
        'Five components primary; three/seven sensitivity and automatic (dim 0) estimation; 24 fits completed. '
        'Every fit uses separate physically masked 4D inputs for its ROI, verified exactly zero outside that ROI. '
        'The joint-BOLD residual subtracts the existing complete fit from centered filtered data; no new regression or temporal filter. '
        'Native FILM residuals have voxelwise prewhitened temporal covariance. '
        'ICA signs/indices are arbitrary; comparisons use one-to-one maximum absolute spatial correlation and retain all correlations/subspace angles. '
        'No task-effect or whole-brain network inference was run.\n\n'+links+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit',type=Path,default=AUDIT)
    parser.add_argument('--jobs',type=int,default=64)
    parser.add_argument('--ica-jobs',type=int,default=18)
    parser.add_argument('--auto-ica-jobs',type=int,default=6,help='Automatic-dimension queue slots; each fit reserves 16 GiB plus shared 8 GiB headroom')
    parser.add_argument('--resume-record',type=Path,help='Reuse the output root of an earlier paired-group-ICA record')
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args()
    core.require(1<=args.jobs<=64 and 1<=args.ica_jobs<=18 and 1<=args.auto_ica_jobs<=6,'Invalid worker limits')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    if args.resume_record:
        previous=json.loads((args.resume_record/'state.json').read_text())
        root=Path(previous['output_root']).resolve()
        core.require(root.is_relative_to((core.ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('paired-group-ica-') and root.is_dir(),'Invalid resume root')
    else:
        root=core.ROOT/'derivatives/insula_ica'/('paired-group-ica-'+stamp);root.mkdir(parents=True)
    record=core.ROOT/'code/insula_ica/records'/('linux2-paired-group-ica-'+stamp);record.mkdir()
    for name in ('errors','figures'):(record/name).mkdir()
    state=dict(state='preflight',runs={},fits={},output_root=str(root),
        started_utc=datetime.now(timezone.utc).isoformat(),code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=core.ROOT,text=True).strip())
    save(record/'state.json',state)
    lock_path=core.ROOT/'code/insula_ica/output/paired-group-ica.lock';lock_path.parent.mkdir(parents=True,exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:run(args,root,record,state)
        except Exception as exc:
            state.update(state='failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        save(record/'state.json',state)
        if args.publish:
            try:publish(record,'Record paired BOLD-scale/native-FILM insula ICA and dimensionality comparisons')
            except Exception as exc:
                (record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
    print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state']=='paired_ICA_complete_review_pending' else 1)


if __name__=='__main__':main()
