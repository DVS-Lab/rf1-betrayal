#!/usr/bin/env python3
"""Functional-cluster mean + separate paired-scale masked ICA; reuse insula fits."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import traceback
import uuid

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code/insula_ica'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import task_dimension_sensitivity as d
import run_network_timecourse_analysis as n
import run_group_ica as g
import resume_group_ica as retry
from run_bold_controls import save, publish

RECORDS=ROOT/'code/functional_cluster_ica/records'
POLICY='functional-127-mean-and-sm0-paired-ICA-task-OLS-v1'
ROI='functional_cluster'
MASK_SHA='66be24139534438a8f2efd49d071afcfa21010df5a578b955917a22b71d234d0'
DIMENSIONS=(3,5,7,0)
SUCCESS='functional_cluster_complete_review_pending'
COLORS=('#68a9cf','#186da2','#e6a658','#bf7340','#88532e')


def read(path):return json.loads(Path(path).read_text())


def atomic_array(path,values):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    with temporary.open('wb') as stream:np.save(stream,values)
    temporary.replace(path)


def preserve(path):
    if path.exists():path.rename(path.with_name(path.name+'.interrupted-'+uuid.uuid4().hex[:6]))


def cache_values(folder,signature,loader):
    """An exact 127-voxel task-preserved cache supplies means and all future ICs."""
    folder.mkdir(parents=True,exist_ok=True);marker=folder/'extraction.json';target=folder/'task_BOLD.npy'
    if marker.is_file():
        old=read(marker);g.core.require(old['signature']==signature and g.core.identity(target)==old['output'],'Functional task cache source/output changed')
        return np.load(target),True
    preserve(target);values=np.asarray(loader(),dtype='float32')
    g.core.require(values.ndim==2 and np.isfinite(values).all() and values.shape[0]==127,'Invalid original functional-mask BOLD cache')
    g.core.require((values.mean(axis=0).std()>1e-8),'Constant functional mean timecourse')
    atomic_array(target,values);save(marker,dict(signature=signature,output=g.core.identity(target)))
    return values,False


def spatial_weights(maps):
    maps=np.asarray(maps,dtype=float)
    g.core.require(maps.ndim==2 and np.isfinite(maps).all(),'Invalid functional ICA maps')
    model=np.column_stack([np.ones(len(maps)),maps-maps.mean(axis=0)])
    g.core.require(np.linalg.matrix_rank(model)==maps.shape[1]+1,'Nonestimable functional spatial basis')
    return np.linalg.pinv(model)[1:]


def fit_response(row,folder,signature,responses,layout,model,cs,labels):
    folder.mkdir(parents=True,exist_ok=True);marker=folder/'complete.json'
    if marker.is_file():
        old=read(marker)
        g.core.require(old['signature']==signature and all(g.core.identity(x['path'])==x for x in old['outputs']),'Completed functional temporal fit changed')
        return old['effects'],True
    preserve(folder/'fit')
    estimates,_=n.fit_timecourses(responses,model,cs,folder/'fit')
    effects=[dict(item,subject=row['subject'],task=row['task'],run=int(row['run']),contrast=label,effect=float(estimates[j,k]))
             for j,label in enumerate(labels) for k,item in enumerate(layout)]
    save(marker,dict(signature=signature,effects=effects,outputs=[g.core.identity(p) for p in (folder/'fit').iterdir() if p.is_file()]))
    return effects,False


def task_context(row,reference,source_output):
    feat=Path(row['active_feat']);design,identities=g.inputs.inspect_source(row,reference)
    settings=g.core.fsl_settings(feat/'design.fsf')
    g.core.require(float(settings['fmri(smooth)'])==0 and settings['fmri(temphp_yn)']=='0' and settings['fmri(ndelete)']=='0','Expected completed unsmoothed/no-filter/no-deletion source')
    bold=n.feat_input(row,settings,reference);saved=g.core.fsl_matrix(feat/'design.con')
    cs,labels=n.contrasts_for(row,settings,design,saved)
    selected=[x for x in identities if Path(x['path']).name in ('design.mat','design.con','confoundevs.txt','mask.nii.gz') or x['path']==str(Path(row['confounds']).resolve())]
    anchors=dict(model_inputs=selected,fsf=g.core.identity(feat/'design.fsf'),FEAT_input=g.core.identity(bold,False),tr=row['tr'])
    old=read(source_output/'runs'/row['key']/'complete.json')['signature']
    g.core.require(all(anchors[k]==old[k] for k in anchors),'Task input/design differs from completed anatomical-ICA task analysis')
    if row['task']=='ugr':design,cs,_=d.pooled.pool_design(design,settings);labels=list(n.BAR_LABELS[:2])
    return design,cs,labels,bold,anchors


def verify_anchors(anchors):
    items=[anchors['fsf'],anchors['FEAT_input'],*anchors['model_inputs']]
    g.core.require(all(g.core.identity(p['path'],p['sha256'] is not None)==p for p in items),'Task sources changed during analysis')


def draw_bars(ax,blocks,names,title):
    for j,(label,color) in enumerate(zip(n.BAR_LABELS,COLORS)):
        rr=[next(r for r in block if r['contrast']==label) for block in blocks]
        y=np.array([float(r['effect']) for r in rr]);lo=np.array([float(r['ci_low']) for r in rr]);hi=np.array([float(r['ci_high']) for r in rr])
        ax.bar(np.arange(len(blocks))+(j-2)*.155,y,width=.145,color=color,label=n.BAR_TITLES[j],yerr=[y-lo,hi-y],capsize=2)
    ax.axhline(0,color='.4',lw=.8);ax.set(xticks=range(len(names)),xticklabels=names,title=title,ylabel='Standardized coefficient / contrast')
    ax.spines[['top','right']].set_visible(False)


def insula_plots(existing,means,folder):
    reference=[r for r in means if r['series']=='functional_mean']
    g.core.require(len(reference)==5,'All five functional-mean bars required')
    for mode in g.inputs.MODES:
        for dim in (3,5,7):
            fig,axes=plt.subplots(3,1,figsize=(15,10),sharey=True)
            for ax,roi in zip(axes,g.core.ROIS):
                parent=[r for r in existing if r['roi']==roi and r['series']=='insula_mean']
                components=[[r for r in existing if r['mode']==mode and r['roi']==roi and int(r['dimensionality'])==dim and int(r['component'])==k] for k in range(1,dim+1)]
                draw_bars(ax,[reference,parent,*components],['Functional\ncluster mean','Insula mean',*['IC'+str(k) for k in range(1,dim+1)]],roi.capitalize()+' anatomical insula')
                ax.axvline(1.5,color='.65',lw=.8,linestyle=':')
            fig.suptitle(f'{mode} | {dim} anatomical-insula components + both mean references',fontsize=14)
            handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.95),ncol=3,frameon=False,fontsize=9)
            fig.text(.5,.025,'N=132; 490 runs. UGR pooled high/low pmods; Trust reciprocate − defect. Pointwise 95% participant CIs. Existing insula estimates reused unchanged.',ha='center',fontsize=8.5)
            fig.tight_layout(rect=(0,.055,1,.89));fig.savefig(folder/f'{mode}_dim-{dim}_insula_ICs_both_means.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,5))
    blocks=[reference]+[[r for r in existing if r['roi']==roi and r['series']=='insula_mean'] for roi in g.core.ROIS]
    draw_bars(ax,blocks,['Functional cluster mean','Left insula mean','Right insula mean','Bilateral insula mean'],'Functional cluster and anatomical insula means')
    ax.legend(ncol=3,loc='upper center',bbox_to_anchor=(.5,-.12),frameon=False,fontsize=8)
    fig.tight_layout(rect=(0,.05,1,1));fig.savefig(folder/'mean_comparison_five_bars.png',dpi=160);plt.close(fig)


def component_plots(groups,means,existing,fits,folder):
    refs=[[r for r in means if r['series']=='functional_mean'],[r for r in existing if r['roi']=='left' and r['series']=='insula_mean']]
    for key,fit in fits.items():
        if fit['state']!='complete':continue
        mode,_,requested=key.split('/');dim=int(requested.removeprefix('dim-'));actual=fit['actual_components']
        for start in range(1,actual+1,8):
            ks=list(range(start,min(start+8,actual+1)))
            blocks=refs+[[r for r in groups if r['mode']==mode and int(r['dimensionality'])==dim and int(r['component'])==k] for k in ks]
            fig,ax=plt.subplots(figsize=(max(12,len(blocks)*1.4),5.8))
            draw_bars(ax,blocks,['Functional mean','Left insula mean',*['IC'+str(k) for k in ks]],f'Functional-cluster ICA | {mode} | requested {dim}; estimated {actual}')
            ax.legend(ncol=3,loc='upper center',bbox_to_anchor=(.5,-.12),frameon=False,fontsize=8)
            fig.text(.5,.012,'Exploratory task-selected ROI. N=132; all runs; pooled UGR + Trust rec−def. Descriptive pointwise 95% CIs, no significance selection.',ha='center',fontsize=8)
            fig.tight_layout(rect=(0,.06,1,1));fig.savefig(folder/f'{mode}_dim-{dim}_functional_task_bars_page-{(start-1)//8+1}.png',dpi=150);plt.close(fig)



def prepare_residuals(row,root,record,original,reference,parent_source):
    """Read retained fit products only; this function has no FILM recovery path."""
    design,before=g.inputs.inspect_source(row,reference)
    parent=read(parent_source/'runs'/(row['key']+'.json'))
    g.core.require(parent['source_identities']==before,'Retained FEAT products differ from completed anatomical ICA input snapshot')
    folder=root/'run_arrays';marker=folder/(row['key']+'_complete.json')
    signature=dict(mask_sha256=MASK_SHA,source_identities=before)
    if marker.is_file():
        old=read(marker)
        g.core.require(old['signature']==signature and all(g.core.identity(v['path'])==v for v in old['arrays']),'Cached functional residual input/output changed')
        save(record/'runs'/(row['key']+'.json'),old['qc']);return old['qc'],True
    for path in folder.glob(row['key']+'_*.npy'):preserve(path)
    native_path=Path(row['native_residual']);g.core.require(native_path.is_file(),'Native residual missing; FILM rerun is not permitted')
    native=nib.load(native_path)
    g.core.require(g.core.same_grid(native,reference) and native.shape==(*reference.shape,len(design)),'Retained native residual geometry changed')
    feat=Path(row['active_feat'])
    y=np.asarray(nib.load(row['filtered_data']).dataobj,dtype='float32')[original].astype(float)
    pe=np.column_stack([np.asarray(nib.load(feat/f'stats/pe{k+1}.nii.gz').dataobj)[original] for k in range(design.shape[1])])
    joint=g.inputs.joint_bold_residual(y,design,pe)
    innovation=np.asarray(native.dataobj,dtype='float32')[original].astype(float)
    fitted=(np.asarray(nib.load(feat/'mask.nii.gz').dataobj)>0)[original]
    sigma=np.asarray(nib.load(feat/'stats/sigmasquareds.nii.gz').dataobj)[original]
    valid=fitted & np.isfinite(sigma) & (sigma>0) & np.isfinite(joint).all(axis=1) & np.isfinite(innovation).all(axis=1)
    valid &= (joint.std(axis=1)>1e-8) & (innovation.std(axis=1)>1e-8)
    g.core.require(valid.any(),'No valid functional residual support')
    paths=[]
    for mode,values in zip(g.inputs.MODES,(joint,innovation)):
        path=folder/(row['key']+'_'+mode+'.npy');atomic_array(path,values.astype('float32'));paths.append(path)
    path=folder/(row['key']+'_valid.npy');atomic_array(path,valid);paths.append(path)
    verify_anchors(dict(fsf=g.core.identity(feat/'design.fsf'),FEAT_input=g.core.identity(row['filtered_data'],False),model_inputs=before))
    g.core.require(before==[g.core.identity(v['path'],v['sha256'] is not None) for v in before],'Retained source changed during residual extraction')
    qc=dict(key=row['key'],valid_voxels=int(valid.sum()),native_recovered_by_film=False,
        source_identities=before,parent_anatomical_ICA_snapshot_verified=True,no_new_regression=True,no_new_temporal_filter=True)
    save(record/'runs'/(row['key']+'.json'),qc);save(marker,dict(signature=signature,qc=qc,arrays=[g.core.identity(p) for p in paths]))
    return qc,False


def build_sessions(rows,root,mask,original,reference):
    slices,affine=g.inputs.crop_geometry(reference,mask);crop=mask[slices]
    g.freeze_nifti(root/'masks'/(ROI+'_crop.nii.gz'),crop.astype('uint8'),affine)
    g.freeze_nifti(root/'reference_crop.nii.gz',reference.get_fdata()[slices].astype('float32'),affine)
    groups=defaultdict(list)
    for row in rows:groups[row['subject'],row['task']].append(row)
    g.core.require(len(groups)==264,'Expected N132 × two task blocks')
    for mode in g.inputs.MODES:
        folder=root/'sessions'/mode/ROI;folder.mkdir(parents=True,exist_ok=True);listing=[];verified=[]
        for (subject,task),rr in sorted(groups.items()):
            rr=sorted(rr,key=lambda r:r['run']);tr=float(rr[0]['tr'])
            g.core.require(all(np.isclose(r['tr'],tr,atol=1e-5) for r in rr),'Mixed within-block TRs; no resampling authorized')
            path=folder/f'sub-{subject}_task-{task}.nii.gz';marker=path.with_name(path.name+'.complete.json')
            signature=dict(run_arrays=[g.core.sha(root/'run_arrays'/(r['key']+'_complete.json')) for r in rr],mask_sha256=g.core.sha(root/'masks'/(ROI+'_crop.nii.gz')),mode=mode,roi=ROI,global_amplitude_multiplier=1000,acquisition_tr_seconds=tr)
            if marker.is_file():
                old=read(marker);check=old['verification'];g.core.require(old['signature']==signature and g.core.identity(path)==check['output'],'Balanced functional 4D cache changed')
            else:
                preserve(path)
                arrays=[np.load(root/'run_arrays'/(r['key']+'_'+mode+'.npy'))[mask[original]] for r in rr]
                values,weights=g.inputs.balanced_session(arrays);data=np.zeros((*crop.shape,values.shape[1]),dtype='float32');data[crop]=values*1000
                g.write_nifti(path,data,affine,tr=tr);check=g.validate_masked_4d(path,nib.load(root/'masks'/(ROI+'_crop.nii.gz')),values.shape[1],tr=tr)
                save(marker,dict(signature=signature,verification=check,within_session_weights=weights))
            listing.append(str(path.resolve()));verified.append(check)
        (folder/'inputs.txt').write_text('\n'.join(listing)+'\n')
        save(folder/'verified-inputs.json',dict(mode=mode,roi=ROI,mask_sha256=g.core.sha(root/'masks'/(ROI+'_crop.nii.gz')),inputs=verified))
        g.verify_ica_inputs(mode,ROI,root)
    return slices


def latest_dimensions():
    candidates=sorted(p for p in (ROOT/'code/insula_ica/records').glob('linux2-task-dimensions-*') if (p/'state.json').is_file())
    g.core.require(bool(candidates),'Completed 3/5/7 task/insula-mean record required')
    source=candidates[-1];state=read(source/'state.json')
    g.core.require(state['state']=='task_dimensions_and_mean_complete_review_pending' and state['policy']==d.POLICY,'Latest anatomical-ICA task analysis is not complete')
    return source


def publish_stage(record,state,enabled,message):
    save(record/'state.json',state)
    if enabled:
        try:
            publish(record,message)
            # A retry can clear a prior publication error after a successful push.
            (record/'publication-error.txt').unlink(missing_ok=True)
        except Exception as exc:
            (record/'publication-error.txt').write_text(str(exc));print('Git publication failed; analysis retained:',exc,flush=True)


def execute(args,record,root,state):
    source=Path(state['source_record']);srcstate=read(source/'state.json');source_output=Path(srcstate['output_root'])
    parent_source=Path(read(source/'map-provenance.json')['ICA_source']);parent_state=read(parent_source/'state.json')
    g.core.require(parent_state['state']=='unsmoothed_ICA_complete_review_pending','Expected completed anatomical ICA source for retained-fit identity checks')
    g.core.require(srcstate['state']=='task_dimensions_and_mean_complete_review_pending' and srcstate['policy']==d.POLICY,'Source task analysis incomplete/changed')
    rows=read(source/'accepted-inputs.json');subjects=sorted({r['subject'] for r in rows})
    cohort=set((ROOT/'code/sublist_n132.txt').read_text().split())
    g.core.require(len(rows)==len({r['key'] for r in rows})==490 and set(subjects)==cohort and len(cohort)==132 and Counter(r['task'] for r in rows)==dict(trust=241,ugr=249),'Manuscript cohort changed')
    g.core.require(all(r['smoothing_mm']==0 for r in rows),'Expected corrected sm-0 source registry')
    reference=nib.load(Path(rows[0]['active_feat'])/'mean_func.nii.gz')
    path=ROOT/'masks/seed-AIns-clusterthresh.nii.gz';g.core.require(g.core.sha(path)==MASK_SHA,'Documented functional mask changed')
    image=nib.load(path);original=image.get_fdata()!=0
    g.core.require(g.core.same_grid(image,reference) and int(original.sum())==127,'Functional seed grid/count changed; no resampling')
    g.freeze_nifti(root/'masks/original_functional_cluster.nii.gz',original.astype('uint8'),reference.affine)
    d.copy_frozen(path,record/'primary_maps/original_functional_cluster.nii.gz')
    shared=dict(policy=POLICY,source_record=str(source),mask=g.core.identity(path),source_groups=g.core.identity(source/'group_task_effects.tsv'),source_registry=g.core.identity(source/'accepted-inputs.json'),
        implementation_sha256={Path(m.__file__).name:g.core.sha(m.__file__) for m in (sys.modules[__name__],n,d,d.pooled,g,g.inputs,retry,g.core)})
    d.freeze_json(record/'provenance.json',shared);d.freeze_json(record/'accepted-inputs.json',rows)
    save(record/'software.json',g.core.software())
    d.copy_frozen(source/'group_task_effects.tsv',record/'reused_insula_group_effects.tsv')
    existing=g.core.read_table(record/'reused_insula_group_effects.tsv');g.core.require(len(existing)==465,'Expected all 465 existing insula group estimates')
    gate=g.inputs.MemoryGate();lock=threading.Lock()
    def update(section,key,result):
        with lock:
            state[section][key]=result;save(record/'state.json',state)
            print(f'{section}: {sum(v["state"]=="complete" for v in state[section].values())}/{490 if section!="fits" else 8} | {key}: {result["state"]}',flush=True)
    def worker(row,section,operation):
        try:
            with gate.reserve(2):value,flags=operation(row)
            result=dict(state='complete',**flags)
            error=record/'errors'/(section+'_'+row['key']+'.txt')
            if error.is_file():error.rename(record/'attempts'/(error.stem+'-'+uuid.uuid4().hex[:6]+'.txt'))
        except Exception as exc:
            value=None;result=dict(state='failed',error=str(exc));(record/'errors'/(section+'_'+row['key']+'.txt')).write_text(traceback.format_exc())
        update(section,row['key'],result);return row['key'],value
    def all_runs(section,operation,pilot=False):
        results={}
        if pilot:
            selected=[rows[0],next(r for r in rows if r['key']=='sub-11031_ses-01_task-ugr_run-1'),next(r for r in rows if r['key']=='sub-10817_ses-01_task-ugr_run-1')]
            results.update(worker(r,section,operation) for r in selected)
            g.core.require(all(v is not None for v in results.values()),section+' pilot failed; no expansion')
            state['pilot_passed']=True;save(record/'state.json',state)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for f in as_completed([pool.submit(worker,r,section,operation) for r in rows if r['key'] not in results]):
                key,value=f.result();results[key]=value
        return results
    def task_values(row):
        model,cs,labels,bold,anchors=task_context(row,reference,source_output)
        signature=dict(anchors,mask_sha256=MASK_SHA)
        values,reused=cache_values(root/'task_values'/row['key'],signature,lambda:np.asarray(nib.load(bold).dataobj,dtype='float32')[original])
        verify_anchors(anchors)
        return values,reused,model,cs,labels,anchors
    def mean(row):
        values,cached,model,cs,labels,anchors=task_values(row)
        layout=[dict(mode='functional_mean',roi=ROI,dimensionality=0,series='functional_mean',component=0)]
        effects,reused=fit_response(row,root/'mean_fits'/row['key'],dict(shared,anchors=anchors,layout=layout),values.astype(float).mean(axis=0)[:,None],layout,model,cs,labels)
        return effects,dict(extraction_reused=cached,fit_reused=reused)
    state['state']='functional_mean_task_fits';save(record/'state.json',state)
    results=all_runs('mean_runs',mean,pilot=True)
    g.core.require(all(v is not None for v in results.values()),'Functional mean fits incomplete; no partial-cohort plot')
    effects=[r for key in sorted(results) for r in results[key]];participant,means=d.aggregate(effects,subjects)
    for name,values in [('functional_mean_run_effects.tsv',effects),('functional_mean_participant_effects.tsv',participant),('functional_mean_group_effects.tsv',means)]:g.core.write_table(record/name,values)
    g.core.write_table(record/'insula_plus_functional_mean_group_effects.tsv',[*existing,*means]);insula_plots(existing,means,record/'figures')
    save(record/'reused-insula-verification.json',dict(estimates_reused=465,new_anatomical_ICA_fits=0,new_anatomical_task_fits=0,original_group_table_sha256=shared['source_groups']['sha256'],copied_group_table_sha256=g.core.sha(record/'reused_insula_group_effects.tsv')))
    state.update(state='functional_mean_plots_ready',mean_plots_ready=True,anatomical_estimates_reused=465)
    publish_stage(record,state,args.publish,'Add functional-cluster mean beside completed insula ICA task bars')
    def residual(row):
        g.core.require(Path(row['native_residual']).is_file(),'Retained native residual missing; no FILM recovery permitted')
        qc,reused=prepare_residuals(row,root,record,original,reference,parent_source)
        return qc,dict(residual_arrays_reused=reused)
    state['state']='extracting_functional_paired_residuals';save(record/'state.json',state)
    results=all_runs('residual_runs',residual,pilot=True)
    g.core.require(all(v is not None for v in results.values()),'Some paired residual extractions failed')
    common=np.ones(127,dtype=bool)
    for row in rows:common &= np.load(root/'run_arrays'/(row['key']+'_valid.npy'))
    mask=np.zeros(reference.shape,dtype=bool);mask[original]=common
    save(record/'coverage.json',dict(original_voxels=127,common_valid_voxels=int(mask.sum()),retained_fraction=float(common.mean()),mean_uses_original_127=True,ICA_uses_common_fitted_finite_nonconstant_support=True,no_subject_or_run_exclusion=True))
    g.core.require(common.mean()>=.8,'Less than 80% common functional-cluster coverage; no automatic exclusions or ICA')
    g.freeze_nifti(root/'masks'/(ROI+'.nii.gz'),mask.astype('uint8'),reference.affine)
    d.copy_frozen(root/'masks'/(ROI+'.nii.gz'),record/'primary_maps/ICA_common_functional_cluster.nii.gz')
    state['state']='building_physically_masked_4d_blocks';save(record/'state.json',state)
    slices=build_sessions(rows,root,mask,original,reference);masks={ROI:mask}
    plans=[(mode,ROI,dim) for mode in g.inputs.MODES for dim in DIMENSIONS]
    save(record/'ica-plan.json',[dict(mode=m,roi=r,dimensionality=z,argv=g.melodic_command(m,r,z,root,variance_normalization=True,preserve_timing=True)) for m,r,z in plans])
    def fit(plan):
        mode,roi,dim=plan;key=f'{mode}/{roi}/dim-{dim}'
        try:
            if dim==0 and Path('/proc/meminfo').is_file():
                total=next(int(line.split()[1])/1024**2 for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemTotal:'))
                g.core.require(total>=40,'Automatic PCA requires 32 GiB plus 8 GiB headroom; fixed fits are retained')
            with gate.reserve(32 if dim==0 else 2):result=retry.retry_fit(plan,root,masks,slices,reference,record,variance_normalization=True,preserve_timing=True,threads=8 if dim==0 else 4)
            d.copy_frozen(Path(result['oriented_maps']),record/'primary_maps'/f'{mode}_dim-{dim}.nii.gz')
            stale=record/'errors'/(key.replace('/','_')+'.txt')
            if stale.is_file():stale.rename(record/'attempts'/(stale.stem+'-'+uuid.uuid4().hex[:6]+'.txt'))
        except Exception as exc:
            result=dict(state='failed',mode=mode,roi=roi,dimensionality=dim,error=str(exc));(record/'errors'/(key.replace('/','_')+'.txt')).write_text(traceback.format_exc())
        update('fits',key,result)
    state['state']='six_fixed_functional_ICA_fits';save(record/'state.json',state)
    with ThreadPoolExecutor(max_workers=args.ica_jobs) as pool:
        for f in as_completed([pool.submit(fit,p) for p in plans if p[2]!=0]):f.result()
    state['state']='two_serial_automatic_functional_ICA_fits';save(record/'state.json',state)
    for plan in plans:
        if plan[2]==0:fit(plan)
    bases={};map_ids={}
    for key,result in state['fits'].items():
        if result['state']!='complete':continue
        data=nib.load(result['oriented_maps']).get_fdata();actual=result['actual_components'];mode=result['mode'];dim=result['dimensionality']
        g.core.require(data.shape==(*reference.shape,actual) and np.count_nonzero(data[~mask])==0,'Functional ICA map support changed')
        bases[key]=(spatial_weights(data[mask]),mode,dim,actual);map_ids[key]=g.core.identity(result['oriented_maps'])
        for start in range(0,actual,12):g.component_figure(data[...,start:start+12],mask,reference.get_fdata(),reference.affine,record/'figures'/f'{mode}_dim-{dim}_spatial_page-{start//12+1}.png',f'Functional cluster | {mode} | requested {dim}',first_component=start+1)
    g.core.require(bool(bases),'No functional ICA completed; functional mean plots remain available')
    comparisons=[]
    for mode in g.inputs.MODES:
        parent=f'{mode}/{ROI}/dim-5'
        if parent not in bases:continue
        primary=nib.load(state['fits'][parent]['oriented_maps']).get_fdata()[mask]
        for dim in (3,7):
            key=f'{mode}/{ROI}/dim-{dim}'
            if key in bases:comparisons.append(dict(mode=mode,sensitivity_dimension=dim,**g.compare_maps(primary,nib.load(state['fits'][key]['oriented_maps']).get_fdata()[mask],labels=('primary','sensitivity'))))
    save(record/'spatial_correspondences.json',comparisons)
    def components(row):
        values,cached,model,cs,labels,anchors=task_values(row)
        retained=values[mask[original]].astype(float)
        effects=[];reuse_flags=[]
        for key,(weights,mode,dim,actual) in bases.items():
            chosen=[dict(mode=mode,roi=ROI,dimensionality=dim,series='component',component=k) for k in range(1,actual+1)]
            result,reused=fit_response(row,root/'component_fits'/mode/f'dim-{dim}'/row['key'],
                dict(shared,anchors=anchors,maps={key:map_ids[key]},layout=chosen),(weights@retained).T,chosen,model,cs,labels)
            effects.extend(result);reuse_flags.append(reused)
        return effects,dict(extraction_reused=cached,fit_reused=all(reuse_flags),basis_fits_reused=sum(reuse_flags),basis_fits_total=len(reuse_flags))
    # Per-basis temporal checkpoints preserve completed task fits when an
    # additional previously failed ICA solution completes on resume.
    save(record/'ICA-task-map-signature.json',map_ids)
    state['state']='functional_component_task_fits';save(record/'state.json',state)
    results=all_runs('component_runs',components,pilot=True)
    g.core.require(all(v is not None for v in results.values()),'Some component task fits failed; no incomplete-cohort summaries')
    effects=[r for key in sorted(results) for r in results[key]];participant,groups=d.aggregate(effects,subjects)
    for name,values in [('functional_IC_run_effects.tsv',effects),('functional_IC_participant_effects.tsv',participant),('functional_IC_group_effects.tsv',groups)]:g.core.write_table(record/name,values)
    component_plots(groups,means,existing,state['fits'],record/'figures')
    state['state']=SUCCESS if all(r['state']=='complete' for r in state['fits'].values()) and len(state['fits'])==8 else 'functional_cluster_partial_ICA_review_errors'
    state.update(ICA_complete=sum(r['state']=='complete' for r in state['fits'].values()),component_task_groups=len(groups),anatomical_estimates_reused=465)
    save(record/'state.json',state)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-record',type=Path);p.add_argument('--resume-record',type=Path);p.add_argument('--new-record',action='store_true')
    p.add_argument('--jobs',type=int,default=32);p.add_argument('--ica-jobs',type=int,default=6);p.add_argument('--publish',action='store_true');args=p.parse_args()
    g.core.require(1<=args.jobs<=64 and 1<=args.ica_jobs<=6,'Expected 1–64 readers and 1–6 fixed ICA slots')
    source=args.source_record or latest_dimensions()
    if not args.resume_record and not args.new_record:
        for candidate in sorted(RECORDS.glob('linux2-functional-cluster-*'),reverse=True):
            if (candidate/'state.json').is_file() and read(candidate/'state.json').get('policy')==POLICY and Path(read(candidate/'state.json')['source_record']).resolve()==source.resolve():args.resume_record=candidate;break
    out=ROOT/'code/insula_ica/output';out.mkdir(parents=True,exist_ok=True)
    with (out/'functional-cluster.lock').open('a') as handle:
        try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:print('Functional-cluster runner already active; no second runner started.');return
        if args.resume_record:
            record=args.resume_record;state=read(record/'state.json');root=Path(state['output_root'])
            g.core.require(state['policy']==POLICY and record.resolve().is_relative_to(RECORDS.resolve()) and root.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('functional-cluster-'),'Invalid resume record/root')
        else:
            stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
            record=RECORDS/('linux2-functional-cluster-'+stamp);root=ROOT/'derivatives/insula_ica'/('functional-cluster-'+stamp)
            record.mkdir(parents=True);root.mkdir(parents=True)
            state=dict(policy=POLICY,state='preflight',source_record=str(source),output_root=str(root),mean_runs={},residual_runs={},fits={},component_runs={},pilot_passed=False)
        for name in ('figures','primary_maps','errors','attempts','runs'):(record/name).mkdir(exist_ok=True)
        for name in ('masks','run_arrays'):(root/name).mkdir(exist_ok=True)
        state.pop('error',None);state.update(jobs=args.jobs,ica_jobs=args.ica_jobs,runner_pid=os.getpid(),last_started_UTC=datetime.now(timezone.utc).isoformat(),code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip());save(record/'state.json',state)
        for path in (record/'failure.txt',):
            if path.is_file():path.rename(record/'attempts'/('failure-'+uuid.uuid4().hex[:6]+'.txt'))
        (record/'report.md').write_text('# Functional cluster: mean and exploratory fractionation\n\nCurrent documented 127-voxel seed-AIns-clusterthresh candidate, not a proven historical PPI mask. Mean: original 127 voxels in exact unresidualized input to completed sm-0 FEAT. New ICA: common valid fitted support of that candidate, minimum 80% retained, no participant/run exclusions. Existing anatomical-insula fits and 465 group estimates remain unchanged.\n\nN132, Trust241/UGR249. Both joint BOLD-scale and native FILM residual representations from retained unsmoothed FEAT fits. Physically masked 4D blocks, per-run SD/run-energy balance, native MELODIC variance normalization enabled, concat, dimensions 3/5/7/automatic. Automatic PCA serialized. No new temporal filtering, task/confound regression, smoothing, fMRIPrep, FEAT or FILM.\n\nNew component timecourses use simultaneous spatial regression inside the functional mask on cached raw FEAT-input BOLD, with a spatial intercept. Responses and full task+nuisance model are z-scored before joint OLS; UGR uses the established pooled high/low fairness slopes, Trust uses C/F/S reciprocate−defect contrasts. Runs averaged within participants, then equal-participant means and descriptive pointwise 95% CIs. Task-selected ROI, exploratory characterization without independent confirmation, permutation tests or significance selection. Component signs set by spatial peak, not task effects; component indices local to each fit.\n\nMeans publish before ICA. All complete raw extractions, temporal fits, residual arrays, masked sessions and ICA fits are checked and reused on resume. New outputs remain under the isolated output_root in state.json; only small records/maps/plots/tables are published. See coverage.json, reused-insula-verification.json, provenance.json, ica-plan.json, fits and errors.\n')
        try:execute(args,record,root,state)
        except Exception as exc:
            state.update(state='functional_cluster_failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        publish_stage(record,state,args.publish,'Publish functional-cluster masked ICA, task bars and reuse checks')
        print('Final state:',state['state'],'| record:',record,flush=True)
        raise SystemExit(0 if state['state']==SUCCESS else 1)


if __name__=='__main__':main()
