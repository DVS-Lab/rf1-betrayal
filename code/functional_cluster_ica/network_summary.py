"""Small-file anatomy and paired task summaries for the fixed functional ICA bases."""
from collections import defaultdict,Counter
from pathlib import Path
import xml.etree.ElementTree as ET
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
from scipy import stats
import run as f

DEFINITIONS={
 'UGR nonsocial':{'ugr_nonsocial':1},'UGR social':{'ugr_social':1},
 'Trust computer':{'trust_C_rec_minus_def':1},'Trust stranger':{'trust_S_rec_minus_def':1},
 'UGR social − nonsocial':{'ugr_social':1,'ugr_nonsocial':-1},
 'Trust stranger − computer':{'trust_S_rec_minus_def':1,'trust_C_rec_minus_def':-1},
 'Trust friend − computer':{'trust_F_rec_minus_def':1,'trust_C_rec_minus_def':-1},
 'Trust social-average − computer':{'trust_F_rec_minus_def':.5,'trust_S_rec_minus_def':.5,'trust_C_rec_minus_def':-1},
 'Trust friend − stranger':{'trust_F_rec_minus_def':1,'trust_S_rec_minus_def':-1},
 'UGR partner average':{'ugr_nonsocial':.5,'ugr_social':.5},
 'Trust partner average':{label:1/3 for label in f.n.BAR_LABELS[2:]}}


def write_table(path,rows):
    """Use LF line endings for new Git-published tables."""
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def load_sources(sm0,sm5):
    sources={'sm0':Path(sm0),'sm5':Path(sm5)};states={k:f.read(p/'state.json') for k,p in sources.items()}
    f.g.core.require(states['sm0']['state']==f.SUCCESS and states['sm5']['state']=='matched_sm5_functional_cluster_complete_review_pending','Both functional ICA arms must have completed')
    rows={arm:f.read(p/'accepted-inputs.json') for arm,p in sources.items()}
    cohort=set((f.ROOT/'code/sublist_n132.txt').read_text().split())
    for rr in rows.values():
        f.g.core.require(len(rr)==len({r['key'] for r in rr})==490 and len(cohort)==132 and {r['subject'] for r in rr}==cohort and Counter(r['task'] for r in rr)==dict(trust=241,ugr=249),'Cohort/run identities differ')
    f.g.core.require([(r['key'],r['tr'],r['nvol']) for r in rows['sm0']]==[(r['key'],r['tr'],r['nvol']) for r in rows['sm5']],'Matched arm timing/order differs')
    masks=[nib.load(p/'primary_maps/ICA_common_functional_cluster.nii.gz') for p in sources.values()]
    mask=masks[0].get_fdata()>0;reference=masks[0]
    f.g.core.require(f.g.core.same_grid(*masks) and np.array_equal(mask,masks[1].get_fdata()>0),'ICA mask differs across smoothing arms')
    original=nib.load(sources['sm0']/'primary_maps/original_functional_cluster.nii.gz').get_fdata()>0
    f.g.core.require(int(original.sum())==127 and int(mask.sum())==117 and np.all(mask<=original),'Frozen functional-mask support changed')
    families={}
    for arm,source in sources.items():
        for mode in f.g.inputs.MODES:
            for dim in (3,5):
                item=states[arm]['fits'][f'{mode}/{f.ROI}/dim-{dim}']
                path=source/'primary_maps'/f'{"sm5_" if arm=="sm5" else ""}{mode}_dim-{dim}.nii.gz'
                f.g.core.require(item['state']=='complete' and f.g.core.sha(path)==item['map_identity']['sha256'],'Published ICA map changed')
                im=nib.load(path);data=im.get_fdata()
                f.g.core.require(data.shape==(*mask.shape,dim) and f.g.core.same_grid(im,reference) and np.isfinite(data).all() and np.count_nonzero(data[~mask])==0,'Invalid functional map support')
                key=f'{arm}_{mode}_dim-{dim}'
                families[key]=dict(key=key,arm=arm,mode=mode,dim=dim,maps=data,map_path=path,map_sha256=f.g.core.sha(path),weights=f.spatial_weights(data[mask]))
    return sources,states,rows,reference,mask,original,families


def paired_summaries(rows,fields):
    """Derive each contrast inside each participant, then summarize participants."""
    grouped=defaultdict(lambda:defaultdict(dict))
    for row in rows:
        key=tuple(row[k] for k in fields);sub=row['subject'];label=row['contrast']
        f.g.core.require(label not in grouped[key][sub],'Duplicate participant task estimate')
        grouped[key][sub][label]=float(row['effect'])
    participants=[];groups=[]
    for key,subjects in sorted(grouped.items()):
        f.g.core.require(len(subjects)==132 and all(set(v)==set(f.n.BAR_LABELS) for v in subjects.values()),'Full paired five-condition cohort required')
        for label,weights in DEFINITIONS.items():
            values=[]
            for subject,estimates in sorted(subjects.items()):
                effect=sum(estimates[name]*weight for name,weight in weights.items());values.append(effect)
                participants.append(dict(zip(fields,key),subject=subject,contrast=label,effect=effect))
            mean=float(np.mean(values));se=float(stats.sem(values));delta=float(stats.t.ppf(.975,131)*se)
            groups.append(dict(zip(fields,key),contrast=label,n=132,effect=mean,se=se,ci_low=mean-delta,ci_high=mean+delta,inference='Descriptive paired-participant contrast; pointwise 95% CI, no multiplicity-adjusted task inference'))
    return participants,groups


def task_summary(sources,record):
    source=sources['sm0'];participants=[]
    mean=f.g.core.read_table(source/'functional_mean_participant_effects.tsv')
    task_source=f.ROOT/'code/insula_ica/records'/Path(f.read(source/'state.json')['source_record']).name
    anatomical=[r for r in f.g.core.read_table(task_source/'participant_task_effects.tsv') if r['series']=='insula_mean' and r['roi']=='left']
    for arm,path in sources.items():
        for r in f.g.core.read_table(path/'functional_IC_participant_effects.tsv'):
            if int(r['dimensionality']) in (3,5):participants.append(dict(r,arm=arm,series='component'))
        for series,values in [('functional_mean',mean),('left_insula_mean',anatomical)]:
            for r in values:participants.append(dict(r,arm=arm,mode='reference',dimensionality='0',series=series,component='0'))
    fields=('arm','mode','dimensionality','series','component');pp,groups=paired_summaries(participants,fields)
    write_table(record/'paired_task_participant_contrasts.tsv',pp);write_table(record/'paired_task_group_contrasts.tsv',groups)
    f.save(record/'task_contrast_definitions.json',dict(definitions=DEFINITIONS,primary_social=['UGR social − nonsocial','Trust stranger − computer'],friend_handling='Separate, with social-average secondary; original five bars retained in source records',orientation='UGR fairness slope and Trust rec-minus-def preserved. Negative social contrasts indicate greater social norm-violation sensitivity.',cross_task_average='No averaging UG and Trust into one coefficient; slopes and binary contrasts differ'))
    for arm in sources:
        for mode in f.g.inputs.MODES:
            for dim in (3,5):
                sets=[[r for r in groups if r['arm']==arm and r['series']==series] for series in ('functional_mean','left_insula_mean')]
                sets += [[r for r in groups if r['arm']==arm and r['mode']==mode and r['dimensionality']==str(dim) and r['component']==str(k)] for k in range(1,dim+1)]
                fig,axes=plt.subplots(2,1,figsize=(14,8))
                for ax,labels,colors,title in [(axes[0],list(DEFINITIONS)[:4],f.COLORS[:2]+(f.COLORS[2],f.COLORS[4]),'Four-condition view: stranger versus computer is matched across tasks'),(axes[1],list(DEFINITIONS)[4:6],('#2489bd','#a65e2e'),'Two primary social-specificity contrasts (computed within participant)')]:
                    count=len(labels);width=.7/count
                    for j,(label,color) in enumerate(zip(labels,colors)):
                        rr=[next(r for r in block if r['contrast']==label) for block in sets];y=np.array([r['effect'] for r in rr]);lo=np.array([r['ci_low'] for r in rr]);hi=np.array([r['ci_high'] for r in rr])
                        ax.bar(np.arange(len(sets))+(j-(count-1)/2)*width,y,width=width*.93,color=color,yerr=[y-lo,hi-y],capsize=2,label=label)
                    ax.axhline(0,color='.4',lw=.8);ax.set(xticks=range(len(sets)),xticklabels=['Functional mean','Left insula mean',*['IC'+str(k) for k in range(1,dim+1)]],ylabel='Standardized coefficient / contrast',title=title)
                    ax.spines[['top','right']].set_visible(False);ax.legend(ncol=2,fontsize=8,frameon=False)
                fig.suptitle(f'{arm} | {mode} | {dim} components',fontsize=14)
                fig.text(.5,.012,'N=132. Original within-task coefficient directions retained; no cross-task difference test. Friend and partner-average summaries remain in tables.\nPointwise paired-participant 95% CIs. Task-selected functional ROI, exploratory. Components are signed spatial patterns, not parcel averages.',ha='center',fontsize=8)
                fig.tight_layout(rect=(0,.065,1,.95),h_pad=1.5);fig.savefig(record/'figures'/f'{arm}_{mode}_dim-{dim}_social_summary.png',dpi=150);plt.close(fig)
    return groups


def anatomy(families,reference,mask,record,fsldir):
    atlaspath=fsldir/'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz';xml=fsldir/'data/atlases/HarvardOxford-Cortical.xml'
    atlas=resample_from_to(nib.load(atlaspath),reference,order=0).get_fdata().astype(int)
    labels={int(x.attrib['index'])+1:x.text for x in ET.parse(xml).findall('./data/label')};labels[0]='Unlabeled'
    bg=resample_from_to(nib.load(fsldir/'data/standard/MNI152_T1_1mm_brain.nii.gz'),reference,order=1).get_fdata();vmax=np.percentile(bg[bg>0],99.5)
    xyz=nib.affines.apply_affine(reference.affine,np.argwhere(mask));coords=[reference.affine[k,3]+np.arange(reference.shape[k])*reference.affine[k,k] for k in range(3)]
    rows=[];detailed=[]
    for key,item in families.items():
        dim=item['dim'];data=item['maps'];fig,axes=plt.subplots(dim,3,figsize=(10,2.5*dim),squeeze=False)
        for k in range(dim):
            values=data[...,k];w=values[mask];energy=w*w;total=energy.sum();peak=np.argwhere(mask)[np.argmax(np.abs(w))];position=xyz[np.argmax(np.abs(w))];center=(xyz*np.abs(w)[:,None]).sum(0)/np.abs(w).sum()
            parts=sorted([dict(label=labels[int(label)],spatial_energy_fraction=float(energy[atlas[mask]==label].sum()/total)) for label in np.unique(atlas[mask])],key=lambda x:-x['spatial_energy_fraction'])
            row=dict(family=key,arm=item['arm'],mode=item['mode'],dimensionality=dim,component=k+1,peak_x=float(position[0]),peak_y=float(position[1]),peak_z=float(position[2]),abs_center_x=float(center[0]),abs_center_y=float(center[1]),abs_center_z=float(center[2]),largest_atlas_energy_label=parts[0]['label'],largest_atlas_energy_fraction=parts[0]['spatial_energy_fraction'])
            rows.append(row);detailed.append(dict(row,atlas_energy=parts))
            for axis,ax in enumerate(axes[k]):
                h,v=((1,2),(0,2),(0,1))[axis];xx,yy=coords[h],coords[v];extent=(xx[0]-reference.header.get_zooms()[h]/2,xx[-1]+reference.header.get_zooms()[h]/2,yy[0]-reference.header.get_zooms()[v]/2,yy[-1]+reference.header.get_zooms()[v]/2);plane=lambda a:np.take(a,int(peak[axis]),axis=axis).T
                ax.imshow(plane(bg),cmap='gray',origin='lower',extent=extent,vmin=0,vmax=vmax,interpolation='nearest')
                ax.imshow(np.ma.masked_where(~plane(mask),plane(values)/np.max(np.abs(w))),cmap='coolwarm',vmin=-1,vmax=1,origin='lower',extent=extent,interpolation='nearest',alpha=.9)
                ax.set_xlim(((8,44),(-63,-21),(-63,-21))[axis]);ax.set_ylim(((-25,23),(-25,23),(8,44))[axis]);ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect('equal');ax.set_title(f'IC{k+1}: {"xyz"[axis]}={position[axis]:.1f} mm',fontsize=9)
            axes[k,0].set_ylabel(parts[0]['label'].replace(', ','\n'),fontsize=8)
        fig.suptitle(key+' | signed component anatomy',fontsize=13)
        fig.text(.5,.009,'All 117 defining voxels displayed; loadings divided by own absolute peak for display (red positive, blue negative).\nNo statistical threshold. Atlas label is largest squared-loading energy, not exclusive parcel membership. MNI coordinates; no map resampling.',ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.045,1,.96));fig.savefig(record/'figures'/f'{key}_anatomy.png',dpi=150);plt.close(fig)
    write_table(record/'component_anatomy.tsv',rows);f.save(record/'component_atlas_loading_profiles.json',detailed)
    f.save(record/'anatomy_sources.json',dict(atlas=f.g.core.identity(atlaspath),labels=f.g.core.identity(xml),interpretation='Squared loading energy across atlas labels; signed weighted patterns, not disjoint parcels'))
    return rows
