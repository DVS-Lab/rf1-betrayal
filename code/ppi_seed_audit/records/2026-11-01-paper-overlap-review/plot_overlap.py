"""Compare original saved whole-brain overlap display across axial slices."""
from pathlib import Path
import hashlib,json,os
import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

record=Path(__file__).resolve().parent;root=record.parents[3]
fsldir=Path(os.environ.get('FSLDIR','/Users/tug87422/fsl'))
bgpath=fsldir/'data/standard/MNI152_T1_1mm_brain.nii.gz'
bgim=nib.as_closest_canonical(nib.load(bgpath));bg=bgim.get_fdata();sources=[];native={};shown={}
for k,name in [('UGR','overlap_figure_ugr_mask.nii.gz'),('Trust','overlap_figure_trust_mask.nii.gz'),('intersection','overlap_figure_intersection_mask.nii.gz'),('seed','seed-AIns-clusterthresh.nii.gz')]:
 path=root/'masks'/name;im=nib.load(path);a=im.get_fdata()!=0;native[k]=a
 sources.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),voxels=int(a.sum())))
 shown[k]=resample_from_to(nib.Nifti1Image(a.astype('uint8'),im.affine),bgim,order=0).get_fdata()>.5
assert np.array_equal(native['UGR']&native['Trust'],native['intersection'])
assert np.array_equal(shown['UGR']&shown['Trust'],shown['intersection'])
# Display categories retain red/blue task-only areas and yellow overlap.
xx=bgim.affine[0,3]+np.arange(bgim.shape[0]);yy=bgim.affine[1,3]+np.arange(bgim.shape[1]);extent=(xx[0]-.5,xx[-1]+.5,yy[0]-.5,yy[-1]+.5)
fig,axes=plt.subplots(3,3,figsize=(12,14));vmax=np.percentile(bg[bg>0],99.5)
cuts=[-16,-13,-10,-7,-4,-1,2,5,8];measure=[]
for ax,z in zip(axes.ravel(),cuts):
 iz=int(round(z-bgim.affine[2,3]));b=bg[:,:,iz].T;u=shown['UGR'][:,:,iz].T;t=shown['Trust'][:,:,iz].T;o=u&t
 rgba=np.zeros((*b.shape,4));rgba[u&~t]=[.9,.16,.13,1];rgba[t&~u]=[.25,.20,.50,1];rgba[o]=[1,.98,0,1]
 ax.imshow(b,cmap='gray',origin='lower',extent=extent,vmin=0,vmax=vmax,interpolation='nearest');ax.imshow(rgba,origin='lower',extent=extent,interpolation='nearest')
 ax.set_xlim(82,-82);ax.set_ylim(-108,84);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);ax.set_title(f'z = {z} mm',fontsize=12)
 ax.text(.035,.04,'R',transform=ax.transAxes,color='white',fontsize=11)
 original_z=int(np.round(nib.affines.apply_affine(np.linalg.inv(im.affine),[0,0,z])[2]))
 measure.append(dict(display_z=z,native_z_index=original_z,native_slice_z=float(nib.affines.apply_affine(im.affine,[0,0,original_z])[2]),whole_brain_overlap_slice_voxels=int(native['intersection'][:,:,original_z].sum()),left_seed_slice_voxels=int(native['seed'][:,:,original_z].sum())))
fig.suptitle('Saved activation masks, shown with the paper’s color convention',fontsize=17,y=.985)
fig.legend(handles=[Patch(color='#e62921',label='UGR only'),Patch(color='#403380',label='Trust only'),Patch(color='#fffa00',label='Trust ∩ UGR')],loc='lower center',bbox_to_anchor=(.5,.033),ncol=3,frameon=False)
fig.text(.5,.014,'Radiological view: participant left is image right. Native masks; nearest-neighbor display, no smoothing or new threshold.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.065,1,.965));fig.savefig(record/'saved_overlap_axial_slices.png',dpi=150);plt.close(fig)
(record/'measurements.json').write_text(json.dumps(dict(sources=sources,slices=measure,note='These are saved overlap_figure masks; no exact correspondence with the published raster figure is assumed.'),indent=2)+'\n')
print(json.dumps(measure,indent=2))

# Display comparison, not an assertion of pixelwise or historical file identity.
import xml.etree.ElementTree as ET
from matplotlib.lines import Line2D
paper=record/'Figure4_original.jpg'
assert paper.is_file(),paper
atlaspath=fsldir/'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz'
labelpath=atlaspath.parent.parent/'HarvardOxford-Cortical.xml'
insid=next(int(e.attrib['index'])+1 for e in ET.parse(labelpath).findall('./data/label') if e.text=='Insular Cortex')
ai=resample_from_to(nib.load(atlaspath),bgim,order=0).get_fdata()==insid
fig,axes=plt.subplots(1,3,figsize=(12,5.5));paperimg=plt.imread(paper)
axes[0].imshow(paperimg);axes[0].set_xlim(130,460);axes[0].set_ylim(562,133);axes[0].set_aspect('equal');axes[0].set_title('Paper: Figure 4A\n(slice coordinate not provided)',fontsize=12)
iz=int(round(-7-bgim.affine[2,3]));b=bg[:,:,iz].T;u=shown['UGR'][:,:,iz].T;t=shown['Trust'][:,:,iz].T
for n,ax in enumerate(axes[1:],1):
 ax.imshow(b,cmap='gray',origin='lower',extent=extent,vmin=0,vmax=vmax,interpolation='nearest')
 rgba=np.zeros((*b.shape,4))
 if n==1:
  rgba[u&~t]=[.9,.16,.13,1];rgba[t&~u]=[.25,.20,.50,1];rgba[u&t]=[1,.98,0,1]
 else:
  s=shown['seed'][:,:,iz].T;rgba[s]=[1,.98,0,1]
  ax.contour(xx,yy,ai[:,:,iz].T,levels=[.5],colors=['#5fe78a'],linewidths=1.3)
 ax.imshow(rgba,origin='lower',extent=extent,interpolation='nearest');ax.set_xlim(82,-82);ax.set_ylim(-108,84);ax.set_aspect('equal');ax.text(.04,.04,'R',transform=ax.transAxes,color='white')
 ax.set_title(('Saved task maps: z = −7 mm\nUGR red · Trust blue · overlap yellow','Selected 127-voxel seed: z = −7 mm\nGreen outline: anatomical insula')[n-1],fontsize=12)
for ax in axes:ax.set_xticks([]);ax.set_yticks([]);[s.set_visible(False) for s in ax.spines.values()]
fig.suptitle('The fuller display closely resembles the paper’s overlap',fontsize=16,y=.98)
fig.text(.5,.032,'Saved-mask native slice: z = −6.66 mm, 19 seed voxels. Earlier centroid slice: z = −3.69 mm, 13 seed voxels.\nReconstruction uses existing masks without changing thresholds. Visual resemblance does not establish exact historical seed identity.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.11,1,.92));fig.savefig(record/'paper_vs_saved_overlap.png',dpi=180);plt.close(fig)
