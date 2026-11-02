"""Read published completed records; summarize all fixed solutions and show aligned dim-5 bars/maps."""
from pathlib import Path
import csv,json,hashlib,sys,os,xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np

record=Path(__file__).resolve().parent;root=record.parents[3]
sys.path.insert(0,str(root/'code/functional_cluster_ica'));import run as f
sources=[record.parent/'linux2-functional-cluster-20261101T200404Z-74f591',record.parent/'linux2-functional-sm5-20261101T211149Z-7aa636']
a,b=sources
load=lambda p:list(csv.DictReader(p.open(),delimiter='\t'))
means=load(a/'functional_mean_group_effects.tsv');old=load(a/'reused_insula_group_effects.tsv');tables=[load(p/'functional_IC_group_effects.tsv') for p in sources]
match=json.loads((b/'smoothing_spatial_correspondence.json').read_text())
mask=nib.load(a/'primary_maps/ICA_common_functional_cluster.nii.gz');inside=mask.get_fdata()>0
original=nib.load(a/'primary_maps/original_functional_cluster.nii.gz').get_fdata()>0
positions=nib.affines.apply_affine(mask.affine,np.argwhere(inside))
fixed=[];scales=[];anatomy=[]
fsldir=Path(os.environ.get('FSLDIR','/Users/tug87422/fsl'));atlasroot=fsldir/'data/atlases';atlaspath=atlasroot/'HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz'
atlas=nib.load(atlaspath);sampled=resample_from_to(atlas,mask,order=0).get_fdata().astype(int)
labels={int(e.attrib['index'])+1:e.text for e in ET.parse(atlasroot/'HarvardOxford-Cortical.xml').findall('./data/label')};labels[0]='Unlabeled'
insula_id=next(k for k,v in labels.items() if v=='Insular Cortex')
images={}
for arm,path in zip(('sm0','sm5'),sources):
 for mode in f.g.inputs.MODES:
  for dim in (3,5,7):
   file=path/'primary_maps'/f'{"sm5_" if arm=="sm5" else ""}{mode}_dim-{dim}.nii.gz';image=nib.load(file);data=image.get_fdata();images[arm,mode,dim]=data
   for k in range(dim):
    w=data[inside,k];weight=np.abs(w);centroid=(weight[:,None]*positions).sum(0)/weight.sum();energy=w*w
    atlas_energy=[dict(label=labels[int(label)],fraction=float(energy[sampled[inside]==label].sum()/energy.sum())) for label in np.unique(sampled[inside])]
    anatomy.append(dict(arm=arm,mode=mode,dimensionality=dim,component=k+1,peak_MNI_mm=positions[np.argmax(np.abs(w))].tolist(),absolute_loading_centroid_MNI_mm=centroid.tolist(),spatial_energy_by_atlas_label=sorted(atlas_energy,key=lambda r:-r['fraction'])))
 for dim in (3,5,7):
  values=f.g.compare_maps(images[arm,'joint_bold',dim][inside],images[arm,'native_film',dim][inside],labels=('bold','film'))
  scales.append(dict(arm=arm,dimensionality=dim,**values))
for item in match:
 vals=np.array([m['absolute_spatial_r'] for m in item['matches']]);fixed.append(dict(mode=item['mode'],dimensionality=item['dimensionality'],matched_components=len(vals),min_abs_r=float(vals.min()),median_abs_r=float(np.median(vals)),max_abs_r=float(vals.max())))
summary=dict(sources=[dict(path=str(p),state=json.loads((p/'state.json').read_text())['state'],group_effects_sha256=hashlib.sha256((p/'functional_IC_group_effects.tsv').read_bytes()).hexdigest()) for p in sources],original_mean_mask_voxels=int(original.sum()),ICA_mask_voxels=int(inside.sum()),whole_ROI_insula_voxels=int(((sampled==insula_id)&original).sum()),ICA_ROI_insula_voxels=int(((sampled==insula_id)&inside).sum()),functional_mean=means,spatial_smoothing_summary=fixed,residual_scale_map_correspondence=scales,component_anatomy=anatomy)
(record/'review-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
# Reference means and every dim-5 IC; sm5 reordered/sign-aligned by maps only.
selected=next(r for r in match if r['mode']=='joint_bold' and r['dimensionality']==5)
ordered=sorted(selected['matches'],key=lambda r:r['sm0_component'])
fig,axes=plt.subplots(3,1,figsize=(13.5,10.7))
f.draw_bars(axes[0],[means,[r for r in old if r['roi']=='left' and r['series']=='insula_mean']],['Functional cluster mean (127 voxels)','Left anatomical insula mean'],'ROI means on identical task-preserved inputs')
f.draw_bars(axes[1],[[r for r in tables[0] if r['mode']=='joint_bold' and r['dimensionality']=='5' and int(r['component'])==k] for k in range(1,6)],['IC'+str(k) for k in range(1,6)],'Unsmoothed functional-cluster ICA: all five components')
blocks=[];names=[]
for m in ordered:
 rr=[dict(r) for r in tables[1] if r['mode']=='joint_bold' and r['dimensionality']=='5' and int(r['component'])==m['sm5_component']]
 if m['sm5_sign_for_sm0_comparison']==-1:
  for r in rr:r['effect']=-float(r['effect']);r['ci_low'],r['ci_high']=-float(r['ci_high']),-float(r['ci_low'])
 blocks.append(rr);names.append(f'5mm IC{m["sm5_component"]}'+(' (flipped)' if m['sm5_sign_for_sm0_comparison']==-1 else '')+f'\n↔ 0mm IC{m["sm0_component"]}; |r|={abs(m["spatial_r"]):.2f}')
f.draw_bars(axes[2],blocks,names,'Matched 5 mm ICA: spatially reordered and sign-aligned, no task-based selection')
for ax in axes[1:]:ax.set_ylim(-.043,.032)
axes[0].set_ylim(-.078,.032)
fig.suptitle('Functional-cluster task responses and smoothing sensitivity',fontsize=16,y=.985)
handles,leg=axes[0].get_legend_handles_labels();fig.legend(handles,leg,loc='upper center',bbox_to_anchor=(.5,.95),ncol=3,frameon=False,fontsize=9)
fig.text(.5,.018,'N=132; 490 runs. Bars: standardized full-model OLS; pointwise 95% participant CIs. No stars or multiplicity-adjusted inference.\nICA uses the same 117 common voxels in both arms. IC responses exclude the spatially uniform signal through the stage-1 spatial intercept.\nMeans are reused raw-input estimates, not a smoothing comparison. Each IC is a signed weighted pattern, not a parcel average.',ha='center',fontsize=8.7)
fig.tight_layout(rect=(0,.085,1,.89),h_pad=2);fig.savefig(record/'five_component_task_summary.png',dpi=170);plt.close(fig)
# Pairwise anatomical display: common sm0-peak cuts, small region, signed maps.
bgpath=fsldir/'data/standard/MNI152_T1_1mm_brain.nii.gz';bgim=nib.as_closest_canonical(nib.load(bgpath));bg=bgim.get_fdata();vmax=np.percentile(bg[bg>0],99.5)
shown_mask=resample_from_to(mask,bgim,order=0).get_fdata()>0
xyz=[bgim.affine[k,3]+np.arange(bgim.shape[k]) for k in range(3)]
fig,axes=plt.subplots(5,6,figsize=(16,12.8));cmap=plt.get_cmap('coolwarm');mappable=None
for i,m in enumerate(ordered):
 raw=images['sm0','joint_bold',5][...,m['sm0_component']-1];other=images['sm5','joint_bold',5][...,m['sm5_component']-1]*m['sm5_sign_for_sm0_comparison']
 peak=positions[np.argmax(np.abs(raw[inside]))];ijk=np.round(nib.affines.apply_affine(np.linalg.inv(bgim.affine),peak)).astype(int)
 for group,values in enumerate((raw,other)):
  values=values/np.abs(values[inside]).max();shown=resample_from_to(nib.Nifti1Image(values,mask.affine),bgim,order=0).get_fdata()
  for axis in range(3):
   ax=axes[i,group*3+axis];h,v=((1,2),(0,2),(0,1))[axis];xx,yy=xyz[h],xyz[v];extent=(xx[0]-.5,xx[-1]+.5,yy[0]-.5,yy[-1]+.5);plane=lambda a:np.take(a,ijk[axis],axis=axis).T
   ax.imshow(plane(bg),cmap='gray',vmin=0,vmax=vmax,extent=extent,origin='lower')
   mappable=ax.imshow(np.ma.masked_where(~plane(shown_mask),plane(shown)),cmap=cmap,vmin=-1,vmax=1,extent=extent,origin='lower',interpolation='nearest',alpha=.9)
   ax.set_xlim(((8,43),(-62,-22),(-62,-22))[axis]);ax.set_ylim(((-24,22),(-24,22),(8,43))[axis]);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
   ax.set_title(f'{"xyz"[axis]}={peak[axis]:.1f}',fontsize=9)
 axes[i,0].set_ylabel(f'0mm IC{m["sm0_component"]}\n↔ 5mm IC{m["sm5_component"]}\nr={abs(m["spatial_r"]):.2f}',fontsize=10)
fig.suptitle('All five functional ICA maps: unsmoothed (left) and matched 5 mm (right)',fontsize=16,y=.985)
fig.text(.5,.012,'Same anatomical cuts within each pair, from the unsmoothed spatial peak. MNI mm; nearest-neighbor display.\nSigned loading divided by its own absolute peak for display (blue negative, red positive). All 117 mask voxels shown; no statistical threshold.\nMaps matched and signs aligned by spatial correlation only. These are overlapping weighted patterns, not five disjoint anatomical parcels.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.065,1,.95),h_pad=.65,w_pad=.7);fig.savefig(record/'five_component_map_comparison.png',dpi=160);plt.close(fig)
print('Smoothing map summaries',json.dumps(fixed[:6]))
for r in anatomy:
 if r['mode']=='joint_bold' and r['dimensionality']==5:print(r['arm'],r['component'],np.round(r['absolute_loading_centroid_MNI_mm'],1).tolist(),r['spatial_energy_by_atlas_label'][:2])
for r in scales:
 if r['dimensionality']==5:print('SCALE agreement',r['arm'],[round(x['absolute_spatial_r'],3) for x in r['matches']])
