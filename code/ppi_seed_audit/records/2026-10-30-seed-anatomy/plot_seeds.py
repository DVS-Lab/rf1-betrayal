"""Plot recovered/current PPI seed candidates; quantify atlas membership on native seed grid."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np


def identity(path):
    return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fsldir',type=Path,required=True);args=p.parse_args()
    record=Path(__file__).resolve().parent;root=record.parents[3]
    old=root/'code/insula_ica/records/2026-10-30-seed-provenance/primary_maps'
    specs=[('Recovered AIns',old/'2aa84f5bd0c7-seed-AIns.nii.gz','#ffae38',142),
           ('Current clusterthresh',root/'masks/seed-AIns-clusterthresh.nii.gz','#36cde2',127),
           ('Recovered AIns-vt',old/'ded74bd4f258-seed-AIns-vt.nii.gz','#f174d4',1)]
    atlasroot=args.fsldir/'data/atlases';ap=atlasroot/'HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz'
    xp=atlasroot/'HarvardOxford-Cortical.xml';atlas=nib.load(ap)
    labels={int(e.attrib['index'])+1:e.text for e in ET.parse(xp).findall('./data/label')};labels[0]='Unlabeled'
    insula=next(k for k,v in labels.items() if v=='Insular Cortex')
    bp=args.fsldir/'data/standard/MNI152_T1_1mm_brain.nii.gz';background=nib.as_closest_canonical(nib.load(bp));bg=background.get_fdata()
    atlas_shown=resample_from_to(atlas,background,order=0).get_fdata()==insula
    rows=[];masks=[];display=[];reference=None
    for name,path,color,count in specs:
        image=nib.load(path);data=image.get_fdata();mask=data!=0
        assert len(image.shape)==3 and np.isfinite(data).all() and int(mask.sum())==count
        if reference is not None:assert image.shape==reference.shape and np.allclose(image.affine,reference.affine)
        reference=image;masks.append(mask)
        if name=='Current clusterthresh':
            assert path.read_bytes()==(old/'4ed77baf3a01-seed-AIns-clusterthresh.nii.gz').read_bytes()
        native=resample_from_to(atlas,image,order=0).get_fdata().astype(int)
        vals,counts=np.unique(native[mask],return_counts=True)
        membership=[dict(label=labels[int(k)],voxels=int(n),percent=float(n/mask.sum()*100)) for k,n in sorted(zip(vals,counts),key=lambda pair:-pair[1])]
        xyz=nib.affines.apply_affine(image.affine,np.argwhere(mask));n=int((native[mask]==insula).sum())
        rows.append(dict(name=name,source=identity(path),voxels=count,voxel_volume_mm3=float(abs(np.linalg.det(image.affine[:3,:3]))),
            centroid_MNI_mm=xyz.mean(0).tolist(),bounds_MNI_mm=[xyz.min(0).tolist(),xyz.max(0).tolist()],
            insula_voxels=n,insula_percent=float(100*n/count),atlas_membership=membership))
        display.append(resample_from_to(nib.Nifti1Image(mask.astype('uint8'),image.affine),background,order=0).get_fdata()>.5)
    overlaps=[]
    for i in range(3):
        for j in range(i+1,3):
            a,b=masks[i],masks[j];shared=int((a&b).sum())
            overlaps.append(dict(first=rows[i]['name'],second=rows[j]['name'],shared_voxels=shared,dice=float(2*shared/(a.sum()+b.sum())),first_only_voxels=int((a&~b).sum()),second_only_voxels=int((b&~a).sum())))
    result=dict(atlas=identity(ap),atlas_labels=identity(xp),background=identity(bp),definition='Harvard-Oxford cortical maximum-probability atlas, threshold 25%; Insular Cortex label, nearest-neighbor sampled at original seed voxel centers. Percent denominator is seed voxel count, not total insula volume.',
        display='Nearest-neighbor seed resampling on 1-mm MNI152 template; no smoothing. Centroid cuts per seed; shared field of view. Original masks are not changed.',seeds=rows,pairwise_overlap=overlaps)
    (record/'anatomy.json').write_text(json.dumps(result,indent=2)+'\n')
    # Common cuts for the two large masks; separate cuts through the vt voxel.
    cuts=[(-44,22,-4),(-44,22,-4),(-34,20,-13)]
    positions=[background.affine[k,3]+np.arange(background.shape[k])*background.affine[k,k] for k in range(3)]
    fig,axes=plt.subplots(3,3,figsize=(12.8,10.8),facecolor='#fafafa')
    vmax=np.percentile(bg[bg>0],99.5)
    for col,((name,path,color,count),shown,coords,row) in enumerate(zip(specs,display,cuts,rows)):
        ijk=np.round(nib.affines.apply_affine(np.linalg.inv(background.affine),coords)).astype(int)
        for axis in range(3):
            ax=axes[axis,col];horizontal,vertical=((1,2),(0,2),(0,1))[axis]
            xx,yy=positions[horizontal],positions[vertical];extent=(xx[0]-.5,xx[-1]+.5,yy[0]-.5,yy[-1]+.5)
            plane=lambda a:np.take(a,ijk[axis],axis=axis).T
            ax.imshow(plane(bg),cmap='gray',origin='lower',extent=extent,vmin=0,vmax=vmax,interpolation='nearest')
            if plane(atlas_shown).any():ax.contour(xx,yy,plane(atlas_shown),levels=[.5],colors=['#6df08e'],linewidths=1.2)
            arr=plane(shown);ax.imshow(np.ma.masked_where(~arr,arr),cmap=ListedColormap([color]),origin='lower',extent=extent,interpolation='nearest',vmin=0,vmax=1,alpha=.92)
            ax.set_xlim(((-38,50),(-72,-8),(-72,-8))[axis]);ax.set_ylim(((-32,32),(-32,32),(-12,52))[axis]);ax.set_aspect('equal')
            ax.set_xlabel(('y (posterior → anterior)','x (left hemisphere)','x (left hemisphere)')[axis],fontsize=9)
            ax.set_ylabel(('z (inferior → superior)','z (inferior → superior)','y (posterior → anterior)')[axis],fontsize=9)
            ax.tick_params(labelsize=8);view=('Sagittal','Coronal','Axial')[axis]
            title=f'{view}: {"xyz"[axis]} = {coords[axis]} mm'
            if axis==0:title=f'{name}\n{count} {"voxel" if count==1 else "voxels"} · {row["insula_voxels"]} in insula ({row["insula_percent"]:.1f}%)\n\n'+title
            ax.set_title(title,fontsize=11,color='#20252b',pad=8)
    fig.suptitle('Where are the three PPI seed candidates?',fontsize=18,y=.985)
    fig.legend(handles=[Patch(facecolor='#ffae38',label='Recovered AIns'),Patch(facecolor='#36cde2',label='Current 127-voxel mask'),Patch(facecolor='#f174d4',label='Recovered AIns-vt'),Line2D([0],[0],color='#49c36b',label='Anatomical insula (HO max-probability, 25%)')],loc='lower center',bbox_to_anchor=(.5,.041),ncol=4,frameon=False,fontsize=9)
    fig.text(.5,.018,'MNI coordinates in mm. Native mask counts; nearest-neighbor display without smoothing. Columns 1–2 share cuts; column 3 cuts through its single voxel.',ha='center',fontsize=8.5,color='#444444')
    fig.tight_layout(rect=(0,.085,1,.96),h_pad=1.1,w_pad=1.7)
    fig.savefig(record/'PPI_seed_anatomy.png',dpi=180);fig.savefig(record/'PPI_seed_anatomy.pdf');plt.close(fig)
    print(json.dumps(dict(seeds=rows,pairwise_overlap=overlaps),indent=2))


if __name__=='__main__':main()
