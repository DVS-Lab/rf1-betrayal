#!/usr/bin/env python3
"""Review corrected UGR group maps against original maps and unchanged Trust."""
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code/insula_ica'))
import insula_ica as core


def touching(mask, seed):
    labels,n=ndimage.label(mask,np.ones((3,3,3)))
    ids=np.unique(labels[seed]);ids=ids[ids!=0]
    return np.isin(labels,ids),len(ids)


def figure(before,after,reference,atlas,background,path):
    seed=before;center=np.round(nib.affines.apply_affine(reference.affine,np.argwhere(seed)).mean(0)).astype(int)
    bg=nib.as_closest_canonical(nib.load(background));data=bg.get_fdata()
    atlas_display=resample_from_to(nib.load(atlas),bg,order=0).get_fdata()==2
    def display(mask):return resample_from_to(nib.Nifti1Image(mask.astype('uint8'),reference.affine),bg,order=0).get_fdata()>.5
    overlays=[(display(before&after),'#f4bc75'),(display(before&~after),'#d94c46'),(display(after&~before),'#2da9e9')]
    vox=np.round(nib.affines.apply_affine(np.linalg.inv(bg.affine),center)).astype(int)
    coords=[bg.affine[k,3]+np.arange(bg.shape[k])*bg.affine[k,k] for k in range(3)]
    fig,axes=plt.subplots(1,3,figsize=(14,5.2))
    for ax,k,label in zip(axes,range(3),('Sagittal','Coronal','Axial')):
        x,y=((1,2),(0,2),(0,1))[k];extent=(coords[x][0]-.5,coords[x][-1]+.5,coords[y][0]-.5,coords[y][-1]+.5)
        plane=lambda a:np.take(a,vox[k],axis=k).T
        ax.imshow(plane(data),origin='lower',extent=extent,cmap='gray',vmin=0,vmax=np.percentile(data[data>0],99.5))
        for mask,color in overlays:
            ax.imshow(np.ma.masked_where(~plane(mask),np.ones(plane(mask).shape)),origin='lower',extent=extent,cmap=ListedColormap([color]),vmin=0,vmax=1,interpolation='nearest')
        if plane(atlas_display).any():ax.contour(coords[x],coords[y],plane(atlas_display),levels=[.5],colors=['#66f28a'],linewidths=1.1)
        ax.set_title(f'{label}: {"xyz"[k]}={center[k]} mm');ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlim(((-100,85),(-85,85),(-85,85))[k]);ax.set_ylim(((-55,90),(-55,90),(-110,85))[k]);ax.set_aspect('equal')
        for position,text in ((.03,'P' if k==0 else 'L'),(.94,'A' if k==0 else 'R')):ax.text(position,.94,text,color='white',transform=ax.transAxes)
    fig.suptitle('Corrected UGR × unchanged Trust: tracked clusterthresh region',fontsize=16,y=.98)
    fig.legend(handles=[Patch(color=c,label=t) for c,t in (('#f4bc75','Shared before/after'),('#d94c46','Original only'),('#2da9e9','Corrected only'))]+[Line2D([0],[0],color='#66f28a',label='Anatomical insula (HO25)')],loc='lower center',bbox_to_anchor=(.5,.08),ncol=4,frameon=False,fontsize=9)
    fig.text(.5,.025,f'Ones model, negative contrasts: left functional-overlap component {before.sum()} → {after.sum()} voxels. Same thresholds; N=132. Nearest-neighbor display only; left is left.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.15,1,.92));fig.savefig(path,dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record',type=Path);parser.add_argument('--atlas',type=Path,required=True);parser.add_argument('--background',type=Path,required=True);args=parser.parse_args()
    record=args.record or sorted((ROOT/'code/ugr_reanalysis/records').glob('linux2-ugr-activation-repair-*'))[-1]
    state=json.loads((record/'state.json').read_text())
    core.require(state['state']=='original_UGR_activation_groups_complete_review_pending' and len(state['cases'])==6 and all(r['state']=='complete' for r in state['cases'].values()),'Expected all six completed repaired groups')
    audit=ROOT/'code/ugr_reanalysis/records'/Path(state['source_audit']).name
    manifest=json.loads((audit/'manifest.json').read_text());by={r['source']:audit/r['copy'] for r in manifest if r.get('copy')};cases=json.loads((audit/'group_cases.json').read_text())
    prefix='/ZPOOL/data/projects/rf1-betrayal';seed_im=nib.load(by[prefix+'/masks/seed-AIns-clusterthresh.nii.gz']);seed=seed_im.get_fdata()!=0
    insula=resample_from_to(nib.load(args.atlas),seed_im,order=0).get_fdata()==2
    comparisons=[];intersections=[];plot_masks=None
    for variant in ('ones','full','covs-vt'):
        trust=next(c for c in cases if c['task']=='trust' and Path(c['path']).parent.name.endswith('-cov-'+variant+'-flame1'))
        for cope in (11,12):
            old=next(c for c in cases if c['task']=='ugr' and Path(c['path']).parent.name.endswith('-cov-'+variant+'-flame1') and f'-cope-{cope}_' in Path(c['path']).name)
            measures=json.loads((record/f'{variant}-cope-{cope}-comparison.json').read_text())
            for sign in (1,2):
                name=f'cluster_mask_zstat{sign}.nii.gz'
                old_im=nib.load(by[old['path']+'/cope1.feat/'+name]);new_im=nib.load(record/'group_maps'/f'{variant}-cope-{cope}-{name}')
                core.require(core.same_grid(old_im,new_im) and core.same_grid(old_im,seed_im),'Different map grids')
                a=old_im.get_fdata()!=0;b=new_im.get_fdata()!=0;measure=next(m for m in measures if m['sign']==sign)
                core.require(int(a.sum())==measure['old_thresholded_voxels'] and int(b.sum())==measure['new_thresholded_voxels'] and int((a&b).sum())==measure['shared_thresholded_voxels'],'Published comparison/count mismatch')
                comparisons.append(dict(variant=variant,cope=cope,**measure,old_anatomical_insula_voxels=int((a&insula).sum()),new_anatomical_insula_voxels=int((b&insula).sum())))
                if cope!=11:continue
                tr=nib.load(by[trust['path']+'/cope1.feat/'+name]);core.require(core.same_grid(tr,seed_im),'Trust grid differs')
                t=tr.get_fdata()!=0;before=a&t;after=b&t;old_region,old_parts=touching(before,seed);new_region,new_parts=touching(after,seed)
                intersections.append(dict(variant=variant,sign=sign,old_intersection_voxels=int(before.sum()),new_intersection_voxels=int(after.sum()),old_insula_voxels=int((before&insula).sum()),new_insula_voxels=int((after&insula).sum()),shared_insula_voxels=int((before&after&insula).sum()),old_seed_region_voxels=int(old_region.sum()),new_seed_region_voxels=int(new_region.sum()),old_region_parts=old_parts,new_region_parts=new_parts,old_seed_retained=int((after&seed).sum())))
                if variant=='ones' and sign==2:plot_masks=(old_region,new_region)
    out=record/'review';out.mkdir(exist_ok=True)
    core.write_table(out/'group_comparisons.tsv',comparisons);core.write_table(out/'task_intersections.tsv',intersections)
    full12=next(c for c in cases if c['task']=='ugr' and '-cov-full-flame1' in Path(c['path']).parent.name and '-cope-12_' in Path(c['path']).name)
    old_mask_im=nib.load(by[full12['path']+'/cope1.feat/cluster_mask_zstat2.nii.gz']);old_mask=old_mask_im.get_fdata()!=0
    xyz=nib.affines.apply_affine(old_mask_im.affine,np.argwhere(old_mask));old_z=nib.load(by[full12['path']+'/cope1.feat/stats/zstat2.nii.gz']).get_fdata()[old_mask];new_z=nib.load(record/'group_maps/full-cope-12-zstat2.nii.gz').get_fdata()[old_mask]
    lost=dict(voxels=int(old_mask.sum()),centroid_mm=xyz.mean(0).tolist(),old_peak_Z=float(old_z.max()),new_peak_Z_in_same_region=float(new_z.max()),old_mean_Z=float(old_z.mean()),new_mean_Z_in_same_region=float(new_z.mean()),new_voxels_above_Z_3p1_in_old_region=int((new_z>=3.1).sum()),note='Old cluster no longer survives same correction; unthresholded effect does not disappear.')
    receipt=dict(source_record=record.name,source_audit=audit.name,source_state=core.identity(record/'state.json'),source_atlas=core.identity(args.atlas),group_comparisons=comparisons,task_intersections=intersections,lost_full_covariate_COPE12_negative_cluster=lost,comparison_scope='Original-model activation only. Same thresholds and unchanged Trust maps; no new seed promotion or PPI result.')
    (out/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    figure(*plot_masks,seed_im,args.atlas,args.background,out/'seed_region_before_after.png')
    print(json.dumps(dict(negative_intersections=[r for r in intersections if r['sign']==2],lost_cluster=lost),indent=2))


if __name__=='__main__':main()
