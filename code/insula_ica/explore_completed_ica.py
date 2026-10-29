#!/usr/bin/env python3
"""Review completed ICA fits now; exploratory residual connectivity and task projections."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
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
from scipy import stats

import run_group_ica as g
from run_bold_controls import save, publish


ROOT=g.core.ROOT
RECORDS=ROOT/'code/insula_ica/records'


def snapshot(source):
    state=json.loads((source/'state.json').read_text());root=Path(state['output_root'])
    g.core.require(root.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('unsmoothed-group-ica-'),'Expected corrected unsmoothed ICA source')
    rows=json.loads((source/'accepted-inputs.json').read_text())
    g.core.require(len(rows)==490 and len({r['subject'] for r in rows})==132,'Expected accepted 490-run/N132 registry')
    completed={}
    for key,item in state.get('fits',{}).items():
        if item.get('state')!='complete' or item['dimensionality']==0:continue
        path=Path(item['oriented_maps'])
        g.core.require(path.resolve().is_relative_to(root.resolve()),'ICA map outside source analysis')
        marker=json.loads((path.parent/'complete.json').read_text())
        g.core.require(marker['map_identity']==g.core.identity(path) and marker['argv']==item['argv'],'Completed ICA map/command changed')
        g.core.require(marker['mask_sha256']==g.core.sha(root/'masks'/(item['roi']+'.nii.gz')) and marker['variance_normalization'],
                       'Completed ICA fixed mask/normalization changed')
        completed[key]=dict(item)
    expected={f'{mode}/{roi}/dim-5' for mode in g.inputs.MODES for roi in g.core.ROIS}
    g.core.require(expected<=set(completed),'All six primary dim-5 fits must be complete; do not select only successful primary effects')
    return root,rows,completed


def weights(maps):
    spatial=np.column_stack([np.ones(len(maps)),maps-maps.mean(axis=0)])
    g.core.require(maps.shape[1]==5 and np.linalg.matrix_rank(spatial)==6,'Nonestimable five-component spatial basis')
    return np.linalg.pinv(spatial)[1:]


def temporal_operator(projector,insula):
    tc=(projector@insula).T;scale=tc.std(axis=0,ddof=1)
    g.core.require(np.isfinite(tc).all() and (scale>1e-8).all(),'Invalid component timecourses')
    standardized=tc/scale;condition=float(np.linalg.cond(standardized))
    g.core.require(np.linalg.matrix_rank(standardized)==5 and condition<1e6,'Ill-conditioned component temporal basis')
    return np.linalg.pinv(standardized),standardized,scale,condition


def joint_coefficients(operator,filtered,design,betas):
    """Exact existing GLS residual projection; no new task/nuisance fit or filter."""
    filtered=np.asarray(filtered,dtype=float);betas=np.asarray(betas,dtype=float)
    return operator@(filtered-filtered.mean(axis=1,keepdims=True)).T-(operator@design)@betas.T


def group_task_tests(run_effects):
    grouped=defaultdict(list)
    for row in run_effects:grouped[(row['mode'],row['roi'],row['subject'],row['task'],row['component'])].append(row['effect'])
    participants=[dict(mode=m,roi=roi,subject=sub,task=task,component=k,runs=len(v),effect=float(np.mean(v)))
                  for (m,roi,sub,task,k),v in sorted(grouped.items())]
    groups=[]
    for mode in g.inputs.MODES:
        for roi in g.core.ROIS:
            for task in g.core.TASKS:
                for k in range(1,6):
                    values=np.array([v['effect'] for v in participants if (v['mode'],v['roi'],v['task'],v['component'])==(mode,roi,task,k)])
                    g.core.require(len(values)==132 and np.isfinite(values).all() and values.std(ddof=1)>0,'Incomplete/nonvarying participant task estimates')
                    mean=float(values.mean());se=float(stats.sem(values));t=mean/se;p=float(2*stats.t.sf(abs(t),131));delta=float(stats.t.ppf(.975,131)*se)
                    groups.append(dict(mode=mode,roi=roi,task=task,component=k,n=132,effect=mean,se=se,ci_low=mean-delta,ci_high=mean+delta,
                        t=t,p=p,p_bonferroni_60=min(1.,60*p),p_bonferroni_30_within_scale=min(1.,30*p),
                        stored_contrast=g.core.TASKS[task][2],direction='As stored; no semantic sign flip',exploratory=True))
    return participants,groups


def review(root,rows,completed,record):
    reference=nib.load(Path(rows[0]['active_feat'])/'mean_func.nii.gz');masks={roi:nib.load(root/'masks'/(roi+'.nii.gz')).get_fdata()>0 for roi in g.core.ROIS}
    anatomy=[];xyz=nib.affines.apply_affine(reference.affine,np.indices(reference.shape).reshape(3,-1).T).reshape(*reference.shape,3)
    (record/'figures').mkdir(exist_ok=True);(record/'primary_maps').mkdir(exist_ok=True)
    for key,item in completed.items():
        roi=item['roi'];mode=item['mode'];dim=item['dimensionality'];im=nib.load(item['oriented_maps']);data=im.get_fdata()
        g.core.require(g.core.same_grid(im,reference) and np.isfinite(data).all() and data.shape[3]==dim,'Unexpected completed map geometry')
        path=record/'figures'/f'{mode}_{roi}_dim-{dim}.png'
        g.component_figure(data,masks[roi],reference.get_fdata(),reference.affine,path,f'Completed unsmoothed ICA | {mode} | {roi} | {dim}')
        if dim!=5:continue
        shutil.copyfile(item['oriented_maps'],record/'primary_maps'/f'{mode}_{roi}_dim-5.nii.gz')
        pos=xyz[masks[roi]]
        for k in range(5):
            v=data[...,k][masks[roi]];a=np.abs(v);center=(pos*a[:,None]).sum(axis=0)/a.sum();peak=pos[np.argmax(a)]
            left=float(np.sum(v[pos[:,0]<0]**2));right=float(np.sum(v[pos[:,0]>0]**2))
            anatomy.append(dict(mode=mode,roi=roi,component=k+1,peak_x=float(peak[0]),peak_y=float(peak[1]),peak_z=float(peak[2]),
                abs_loading_center_x=float(center[0]),abs_loading_center_y=float(center[1]),abs_loading_center_z=float(center[2]),
                right_minus_left_energy_ratio=(right-left)/(right+left),
                note='Descriptive signed component loadings; not disjoint parcels. Largest absolute peak oriented positive.'))
    g.core.write_table(record/'component_anatomy.tsv',anatomy)
    comparisons=[]
    for mode in g.inputs.MODES:
        for roi in g.core.ROIS:
            main=nib.load(completed[f'{mode}/{roi}/dim-5']['oriented_maps']).get_fdata()[masks[roi]]
            for dim in (3,7):
                key=f'{mode}/{roi}/dim-{dim}'
                if key in completed:
                    other=nib.load(completed[key]['oriented_maps']).get_fdata()[masks[roi]]
                    comparisons.append(dict(mode=mode,roi=roi,sensitivity_dimension=dim,
                        **g.compare_maps(main,other,labels=('primary','sensitivity'))))
    save(record/'fixed_dimension_comparisons.json',comparisons)
    primary={}
    for mode in g.inputs.MODES:
        for roi in g.core.ROIS:
            item=completed[f'{mode}/{roi}/dim-5'];maps=nib.load(item['oriented_maps']).get_fdata()
            primary[mode,roi]=dict(item=item,mask=masks[roi],maps=maps,weights=weights(maps[masks[roi]]))
    save(record/'between_scale_primary_comparisons.json',[dict(roi=roi,**g.compare_maps(
        primary['joint_bold',roi]['maps'][masks[roi]],primary['native_film',roi]['maps'][masks[roi]])) for roi in g.core.ROIS])
    return reference,masks,primary


def task_figures(groups,record):
    for mode in g.inputs.MODES:
        fig,axes=plt.subplots(3,2,figsize=(11,10))
        for i,roi in enumerate(g.core.ROIS):
            for j,task in enumerate(g.core.TASKS):
                values=[r for r in groups if (r['mode'],r['roi'],r['task'])==(mode,roi,task)];ax=axes[i,j]
                y=np.array([r['effect'] for r in values]);lo=np.array([r['ci_low'] for r in values]);hi=np.array([r['ci_high'] for r in values])
                ax.errorbar(np.arange(1,6),y,yerr=[y-lo,hi-y],fmt='o',capsize=3);ax.axhline(0,color='.5',lw=1)
                for k,r in enumerate(values):
                    if r['p_bonferroni_60']<.05:ax.annotate('*',(k+1,y[k]),xytext=(5,5),textcoords='offset points')
                ax.set(title=f'{roi} | {task}: {g.core.TASKS[task][2]}',xlabel='Component (within this ROI/scale)',ylabel='Projected stored COPE / spatial-map unit',xticks=range(1,6))
        fig.suptitle(f'{mode}: exploratory task-contrast projections, N=132\n95% pointwise participant CIs; * Bonferroni across 60 tests; signs as stored',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.94));fig.savefig(record/'figures'/f'{mode}_task_projections.png',dpi=130);plt.close(fig)


def connectivity_figure(primary,mean,roi,mode,wb,reference,path):
    bg=reference.get_fdata();mask=primary['mask'];maps=primary['maps'];outside=wb&~mask
    fig,axes=plt.subplots(5,2,figsize=(11,15),facecolor='black')
    for k in range(5):
        for col,(data,valid) in enumerate(((maps[...,k],mask),(mean[...,k],outside))):
            vals=np.abs(data[valid]);threshold=float(np.percentile(vals,70));scale=max(float(np.percentile(vals,98)),1e-8)
            peak=np.argmax(np.where(valid,np.abs(data),-1));z=np.unravel_index(peak,data.shape)[2]
            plane=np.ma.masked_where(~valid[:,:,z].T | (np.abs(data[:,:,z].T)<=threshold),data[:,:,z].T)
            ax=axes[k,col];ax.imshow(bg[:,:,z].T,origin='lower',cmap='gray',vmin=0,vmax=np.percentile(bg[bg>0],99))
            ax.imshow(plane,origin='lower',cmap='coolwarm',vmin=-scale,vmax=scale,alpha=.85);ax.axis('off')
            ax.set_title(f'IC {k+1}: '+('insula loading' if col==0 else 'outside-ROI background coefficient')+f'; z={reference.affine[2,3]+z*reference.affine[2,2]:g} mm',color='white',fontsize=9)
    fig.suptitle(f'{mode} | {roi}: component anatomy and descriptive background connectivity\nLargest 30% |values| displayed; separate scales; no voxelwise significance threshold',color='white',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(path,dpi=115,facecolor='black');plt.close(fig)


def run(args,record,output,state):
    root,rows,completed=snapshot(args.source_record)
    save(record/'source_snapshot.json',dict(source_record=str(args.source_record),analysis_root=str(root),fits=completed,
        selection='All completed fixed fits for review; all six dim-5 bases fixed before examining task/connectivity results.'))
    reference,masks,primary=review(root,rows,completed,record)
    policy=dict(source_record=str(args.source_record),analysis_scope='Exploratory completed primary ICA maps; automatic fits not required.',
        wholebrain_connectivity='Same residual representation in spatial and temporal regression; all five component timecourses entered simultaneously.',
        extra_task_or_nuisance_regression=False,extra_temporal_filter=False,new_ICA=False,new_FEAT=False,
        task_projection='Same frozen multivariate spatial weights applied to the current corrected sm-0 FEAT COPE.',
        task_directions='As stored: Trust rec-def; UGR offer pmod. No inferred unfairness sign flip.',
        multiplicity='60 tests across both residual bases, ROIs, components and tasks; within-scale 30-test values also saved.',
        participant_aggregation='Equal run means within participant/task, then equal participant means. Connectivity pooling gives each task equal weight.',
        task_inference='Exploratory; same-cohort map discovery, no independent validation or task interaction claimed.',
        native_FIlM='Innovation-scale background coefficients, distinct from original-BOLD covariance.',
        seed_ROI_excluded_from_connectivity_display=True,network_names='Anatomical atlas target descriptions only; no reference-network assignment.')
    save(record/'policy.json',policy)
    state.update(state='completed_fixed_map_review_published',source_fits=len(completed));save(record/'state.json',state)
    (record/'README.md').write_text('# Completed unsmoothed ICA exploration\n\n'
        'Fixed-fit spatial figures and primary component maps are available now. Automatic ICA continues independently. '
        'Whole-brain connectivity and task projections are pending unless state.json reports their completion. '
        'No network/task findings are inferred from convergence alone. See policy.json and component_anatomy.tsv.\n\n'+
        '\n'.join(f'- [{v["mode"]}, {v["roi"]}, {v["dimensionality"]} components](figures/{v["mode"]}_{v["roi"]}_dim-{v["dimensionality"]}.png)' for v in completed.values())+'\n')
    if args.publish:
        try:publish(record,'Publish completed unsmoothed ICA component review while automatic fits continue')
        except Exception as exc:(record/'review-publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
    if args.review_only:return
    state['state']='projecting_corrected_task_contrasts';save(record/'state.json',state)
    effects=[];cope_identities=[]
    for row in rows:
        path=Path(row['active_feat'])/f'stats/cope{row["contrast_index"]}.nii.gz';identity=g.core.identity(path,False)
        im=nib.load(path);g.core.require(g.core.same_grid(im,reference),'Task COPE grid differs')
        data=im.get_fdata()
        for (mode,roi),current in primary.items():
            g.core.require(np.isfinite(data[current['mask']]).all(),'Nonfinite task COPE within insula')
            estimate=current['weights']@data[current['mask']]
            for k,value in enumerate(estimate):effects.append(dict(mode=mode,roi=roi,subject=row['subject'],task=row['task'],run=row['run'],component=k+1,
                effect=float(value),source_cope=str(path),contrast_index=row['contrast_index']))
        g.core.require(g.core.identity(path,False)==identity,'Task COPE changed during projection');cope_identities.append(identity)
    save(record/'task_cope_provenance.json',cope_identities)
    participant,groups=group_task_tests(effects)
    g.core.write_table(record/'run_task_projections.tsv',effects);g.core.write_table(record/'participant_task_projections.tsv',participant)
    g.core.write_table(record/'group_task_projections.tsv',groups);task_figures(groups,record)
    state.update(state='task_projection_review_ready_connectivity_pending',tests=60,
        task_tests_p_bonferroni_60_below_05=sum(v['p_bonferroni_60']<.05 for v in groups));save(record/'state.json',state)
    with (record/'README.md').open('a') as stream:
        stream.write('\n\n## Task-projection review ready\n\n'
            f'All 60 planned projections are available; {state["task_tests_p_bonferroni_60_below_05"]} have Bonferroni-60 p<0.05. '
            'Same-cohort exploratory estimates, not independent confirmation. Direction is as stored in each task contrast. '
            'Connectivity remains pending at this checkpoint.\n\n'
            '- [Joint-BOLD task projections](figures/joint_bold_task_projections.png)\n'
            '- [Native-FILM task projections](figures/native_film_task_projections.png)\n'
            '- [All task estimates and CIs](group_task_projections.tsv)\n')
    if args.publish:
        try:publish(record,'Publish all primary ICA task-contrast projections before whole-brain connectivity finishes')
        except Exception as exc:(record/'task-publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
    g.core.require(shutil.disk_usage(output).free>32*1024**3,'Need 32 GiB free for retained run/participant connectivity maps')
    fsldir=Path(os.environ['FSLDIR']);atlaspath=fsldir/'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'
    xml=fsldir/'data/atlases/HarvardOxford-Cortical.xml'
    atlas_source=json.loads((args.source_record/'atlas.json').read_text())
    g.core.require(g.core.sha(atlaspath)==atlas_source['atlas']['sha256'] and g.core.sha(xml)==atlas_source['xml']['sha256'],
                   'Atlas changed since residual array construction; anatomical voxel ordering must remain identical')
    anatomical,_,_=g.core.anatomical_masks(atlaspath,xml,reference)
    union=anatomical['bilateral'];wb=np.ones(reference.shape,bool);fitted_union=np.zeros(reference.shape,bool)
    for row in rows:
        feat=Path(row['active_feat']);mask=nib.load(feat/'mask.nii.gz');sigma=nib.load(feat/'stats/sigmasquareds.nii.gz')
        g.core.require(g.core.same_grid(mask,reference) and g.core.same_grid(sigma,reference),'Whole-brain source grid differs')
        values=sigma.get_fdata();valid=(mask.get_fdata()>0)&np.isfinite(values)&(values>0);wb&=valid;fitted_union|=valid
    g.core.require(wb.sum()>1000,'Insufficient common whole-brain support')
    for roi in g.core.ROIS:g.core.require((masks[roi]&~wb).sum()==0,'ICA mask extends beyond common fitted whole-brain support')
    g.freeze_nifti(output/'wholebrain_common_mask.nii.gz',wb.astype('uint8'),reference.affine)
    (record/'group_maps').mkdir(exist_ok=True)
    g.core.save_image(record/'group_maps/wholebrain_common_mask.nii.gz',wb.astype('uint8'),reference)
    save(record/'wholebrain_coverage.json',dict(common_voxels=int(wb.sum()),union_voxels=int(fitted_union.sum()),
        common_fraction_of_union=float(wb.sum()/fitted_union.sum()),
        definition='Intersection of positive finite fitted sigma and FEAT masks over all 490 runs. No resampling of data. Outside-mask zeros are missing coverage, not null connectivity.'))
    for name in ('runs','participants','group'):(output/name).mkdir(exist_ok=True)
    gate=g.inputs.MemoryGate();lock=threading.Lock()
    source_signature=dict(maps={f'{m}/{r}':primary[m,r]['item']['map_identity'] for m,r in primary},common_mask=g.core.identity(output/'wholebrain_common_mask.nii.gz'))
    state['state']='running_background_connectivity_and_task_projections';save(record/'state.json',state)
    def worker(row):
        key=row['key'];folder=output/'runs'/key;folder.mkdir(exist_ok=True);marker=folder/'complete.json'
        try:
            source_marker=root/'run_arrays'/(key+'_complete.json')
            checkpoint=json.loads(source_marker.read_text())
            g.core.require(all(g.core.identity(x['path'])==x for x in checkpoint['arrays']),'Frozen insula residual arrays changed')
            g.core.require(all(g.core.identity(x['path'],x['sha256'] is not None)==x for x in checkpoint['qc']['source_identities']),
                           'Frozen fitted model/data inputs changed')
            signature=dict(source_signature,run=g.core.sha(source_marker),
                model=g.core.identity(Path(row['active_feat'])/'design.mat'),
                cope=g.core.identity(Path(row['active_feat'])/f'stats/cope{row["contrast_index"]}.nii.gz',False))
            if marker.is_file():
                old=json.loads(marker.read_text());g.core.require(old['signature']==signature and
                    all(g.core.identity(x['path'])==x for x in old['outputs']),'Completed exploration inputs/outputs changed')
                result=old['result']
            else:
                g.core.require(not any(folder.iterdir()),'Partial exploration run exists; preserve and inspect: '+str(folder))
                with gate.reserve(4):
                    result=analyze_run(row,root,folder,primary,union,wb,reference)
                paths=list(folder.glob('*connectivity.nii.gz'))+list(folder.glob('*timecourses.tsv'))
                save(marker,dict(signature=signature,result=result,outputs=[g.core.identity(p) for p in paths]))
            status=dict(state='complete')
        except Exception as exc:
            result=None;status=dict(state='failed',error=str(exc));(record/'errors'/(key+'.txt')).write_text(traceback.format_exc())
        with lock:
            state['runs'][key]=status;save(record/'state.json',state)
            if len(state['runs'])%25==0 or status['state']=='failed':print(f'Exploration runs {len(state["runs"])}/490 | {key}: {status["state"]}',flush=True)
        return result
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results=[future.result() for future in as_completed([pool.submit(worker,row) for row in rows])]
    g.core.require(all(r is not None for r in results),'Some exploratory runs failed; no incomplete-cohort tests summarized')
    for identity in cope_identities:g.core.require(g.core.identity(identity['path'],False)==identity,'Task COPE changed during connectivity stage')
    state['state']='aggregating_participant_connectivity';save(record/'state.json',state)
    aggregate(rows,output,record,primary,wb,reference,atlaspath,xml)
    state.update(state='exploratory_connectivity_and_task_review_ready',tests=60,
                 task_tests_p_bonferroni_60_below_05=sum(v['p_bonferroni_60']<.05 for v in groups))
    save(record/'state.json',state)
    with (record/'README.md').open('a') as stream:
        stream.write('\n\nConnectivity and all 60 task projections completed for N=132 / 490 runs. '
            f'{state["task_tests_p_bonferroni_60_below_05"]}/60 task tests have p Bonferroni-60 <0.05. '
            'All estimates, CIs and null results are in group_task_projections.tsv. '
            'Directions use the stored contrasts; no unfairness-direction sign was inferred. '
            'Connectivity images are descriptive means, not thresholded group inference. '
            'Atlas regional targets are anatomical summaries, not independent reference-network identification. '
            'Native-FILM maps describe innovation-scale coupling; joint-BOLD maps describe original-scale residual coupling. '
            'These are same-cohort exploratory associations; no unique betrayal encoding, causal role, independent validation, or formal between-task/hemisphere effect is established.\n\n'+
            '\n'.join(f'- [{mode}, {roi}: anatomy and outside-ROI connectivity](figures/{mode}_{roi}_anatomy_connectivity.png)' for mode in g.inputs.MODES for roi in g.core.ROIS)+'\n')


def analyze_run(row,source,folder,primary,union,wb,reference):
    feat=Path(row['active_feat']);design=g.core.fsl_matrix(feat/'design.mat');cope=nib.load(feat/f'stats/cope{row["contrast_index"]}.nii.gz').get_fdata()
    effects=[];diagnostics=[]
    for mode in g.inputs.MODES:
        operators={}
        insula=np.load(source/'run_arrays'/(row['key']+'_'+mode+'.npy'))
        for roi in g.core.ROIS:
            current=primary[mode,roi];operator,tc,scale,condition=temporal_operator(current['weights'],insula[current['mask'][union]])
            operators[roi]=operator;np.savetxt(folder/f'{mode}_{roi}_timecourses.tsv',tc,delimiter='\t')
            estimate=current['weights']@cope[current['mask']]
            for k,value in enumerate(estimate):effects.append(dict(mode=mode,roi=roi,subject=row['subject'],task=row['task'],run=row['run'],
                component=k+1,effect=float(value),source_cope=str(feat/f'stats/cope{row["contrast_index"]}.nii.gz'),contrast_index=row['contrast_index']))
            diagnostics.append(dict(mode=mode,roi=roi,temporal_condition=condition,timecourse_SD=scale.tolist()))
        matrices={roi:np.empty((int(wb.sum()),5),dtype='float32') for roi in g.core.ROIS}
        if mode=='native_film':
            data=np.asarray(nib.load(row['native_residual']).dataobj,dtype='float32')[wb]
            g.core.require(data.shape[1]==row['nvol'] and np.isfinite(data).all(),'Invalid whole-brain native residuals')
            for roi,operator in operators.items():matrices[roi][:]=(operator@data.T).T
            del data
        else:
            data=np.asarray(nib.load(row['filtered_data']).dataobj,dtype='float32')[wb]
            betas=np.column_stack([np.asarray(nib.load(feat/f'stats/pe{k+1}.nii.gz').dataobj,dtype='float32')[wb] for k in range(design.shape[1])])
            g.core.require(np.isfinite(data).all() and np.isfinite(betas).all(),'Nonfinite whole-brain joint-fit inputs')
            for start in range(0,len(data),4096):
                stop=start+4096
                for roi,operator in operators.items():matrices[roi][start:stop]=joint_coefficients(operator,data[start:stop],design,betas[start:stop]).T
            del data,betas
        for roi,matrix in matrices.items():
            g.core.require(np.isfinite(matrix).all(),'Nonfinite background coefficients')
            image=np.zeros((*reference.shape,5),dtype='float32');image[wb]=matrix
            g.core.save_image(folder/f'{mode}_{roi}_connectivity.nii.gz',image,reference)
    return dict(effects=effects,diagnostics=diagnostics)


def aggregate(rows,output,record,primary,wb,reference,atlaspath,xml):
    labels={int(x.attrib['index'])+1:(x.text or '').strip() for x in ET.parse(xml).findall('.//label')}
    atlas=resample_from_to(nib.load(atlaspath),(reference.shape,reference.affine),order=0).get_fdata().astype(int)
    regional=[];subjects=sorted({r['subject'] for r in rows})
    for mode in g.inputs.MODES:
        for roi in g.core.ROIS:
            total=np.zeros((*reference.shape,5),dtype='float64');task_totals={t:np.zeros_like(total) for t in g.core.TASKS}
            for sub in subjects:
                task_means=[]
                for task in g.core.TASKS:
                    rr=[r for r in rows if r['subject']==sub and r['task']==task];g.core.require(bool(rr),'Missing participant/task')
                    mean=sum(nib.load(output/'runs'/r['key']/f'{mode}_{roi}_connectivity.nii.gz').get_fdata(dtype='float32')/len(rr) for r in rr)
                    g.core.save_image(output/'participants'/f'{mode}_{roi}_sub-{sub}_task-{task}.nii.gz',mean.astype('float32'),reference)
                    task_means.append(mean);task_totals[task]+=mean
                total+=np.mean(task_means,axis=0)
            mean=(total/len(subjects)).astype('float32');g.core.save_image(output/'group'/f'{mode}_{roi}_pooled_mean.nii.gz',mean,reference)
            # Group spatial maps are small final artifacts; no subject 4D series
            # or run maps are copied into Git records.
            (record/'group_maps').mkdir(exist_ok=True)
            g.core.save_image(record/'group_maps'/f'{mode}_{roi}_pooled_mean.nii.gz',mean,reference)
            for task,values in task_totals.items():g.core.save_image(output/'group'/f'{mode}_{roi}_{task}_mean.nii.gz',(values/len(subjects)).astype('float32'),reference)
            connectivity_figure(primary[mode,roi],mean,roi,mode,wb,reference,record/'figures'/f'{mode}_{roi}_anatomy_connectivity.png')
            outside=wb&~primary[mode,roi]['mask']
            for index,label in labels.items():
                mask=(atlas==index)&outside
                if mask.sum()<20:continue
                for k,value in enumerate(mean[mask].mean(axis=0)):
                    regional.append(dict(mode=mode,roi=roi,component=k+1,target=label,outside_ROI_voxels=int(mask.sum()),
                        mean_background_coefficient=float(value),descriptive_only=True))
    g.core.write_table(record/'anatomical_connectivity_targets.tsv',regional)
    ranks=[]
    for mode in g.inputs.MODES:
        for roi in g.core.ROIS:
            for k in range(1,6):
                values=[r for r in regional if (r['mode'],r['roi'],r['component'])==(mode,roi,k)]
                ranks.append(dict(mode=mode,roi=roi,component=k,
                    largest_positive=[r for r in sorted(values,key=lambda r:-r['mean_background_coefficient']) if r['mean_background_coefficient']>0][:5],
                    largest_negative=[r for r in sorted(values,key=lambda r:r['mean_background_coefficient']) if r['mean_background_coefficient']<0][:5],
                    note='Descriptive anatomical target means outside this seed ROI; no reference-network or group-significance claim.'))
    save(record/'ranked_anatomical_targets.json',ranks)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-record',type=Path)
    parser.add_argument('--resume-record',type=Path);parser.add_argument('--review-only',action='store_true');parser.add_argument('--jobs',type=int,default=8)
    parser.add_argument('--publish',action='store_true');args=parser.parse_args();g.core.require(1<=args.jobs<=8,'Use 1–8 concurrent whole-brain readers while automatic PCA runs')
    if args.resume_record:
        prior=json.loads((args.resume_record/'state.json').read_text());args.source_record=Path(prior['source_record']);output=Path(prior['output_root'])
    else:output=None
    if args.source_record is None:
        choices=sorted(p for p in RECORDS.glob('linux2-unsmoothed-group-ica-*') if (p/'state.json').is_file())
        g.core.require(bool(choices),'No current unsmoothed ICA record found');args.source_record=choices[-1]
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    if output is None:output=ROOT/'derivatives/insula_ica'/('completed-ICA-exploration-'+stamp)
    g.core.require(output.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve()) and output.name.startswith('completed-ICA-exploration-'),'Invalid exploration output root')
    lockdir=ROOT/'code/insula_ica/output';lockdir.mkdir(exist_ok=True)
    with (lockdir/'completed-ICA-exploration.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        record=RECORDS/('linux2-completed-ICA-exploration-'+stamp);record.mkdir();(record/'errors').mkdir();output.mkdir(parents=True,exist_ok=True)
        state=dict(state='snapshotting_completed_fits',source_record=str(args.source_record),output_root=str(output),runs={})
        save(record/'state.json',state)
        try:run(args,record,output,state)
        except Exception as exc:
            state.update(state='failed',error=str(exc));save(record/'state.json',state)
            (record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        if args.publish:
            try:publish(record,'Publish exploratory primary-component task projections and whole-brain residual connectivity')
            except Exception as exc:(record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
        print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state'] in ('completed_fixed_map_review_published','exploratory_connectivity_and_task_review_ready') else 1)


if __name__=='__main__':main()
