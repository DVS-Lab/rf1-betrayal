"""Display the saved functional seed and verify its two-task intersection support."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
from scipy import ndimage


def identity(path):
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--background', type=Path, required=True)
    args = parser.parse_args()
    record = Path(__file__).resolve().parent; root = record.parents[3]
    names = ('seed-AIns-clusterthresh.nii.gz', 'ugr-trust_intersection.nii.gz',
             'overlap_figure_ugr_mask.nii.gz', 'overlap_figure_trust_mask.nii.gz')
    images = {name:nib.load(root/'masks'/name) for name in names}
    seed_image = images[names[0]]
    assert all(im.shape == seed_image.shape and np.allclose(im.affine, seed_image.affine) for im in images.values())
    masks = {name:im.get_fdata() != 0 for name, im in images.items()}
    seed, intersection, ugr, trust = (masks[name] for name in names)
    assert np.array_equal(intersection, ugr & trust)
    labels, count = ndimage.label(intersection, np.ones((3, 3, 3)))
    seed_labels = np.unique(labels[seed]); assert len(seed_labels) == 1 and seed_labels[0] != 0
    assert np.array_equal(seed, labels == seed_labels[0])
    xyz = nib.affines.apply_affine(seed_image.affine, np.argwhere(seed))
    voxel_volume = float(abs(np.linalg.det(seed_image.affine[:3, :3])))
    atlas_root = args.background.parent.parent/'atlases'
    atlas_path = atlas_root/'HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz'
    label_path = atlas_root/'HarvardOxford-Cortical.xml'
    atlas = nib.load(atlas_path)
    atlas_labels = {int(item.attrib['index'])+1:item.text for item in ET.parse(label_path).getroot().findall('./data/label')}
    insula_index = next(k for k,v in atlas_labels.items() if v == 'Insular Cortex')
    atlas_labels[0] = 'Unlabeled'
    native_labels = resample_from_to(atlas, seed_image, order=0).get_fdata().astype(int)
    counts = np.bincount(native_labels[seed])
    atlas_overlap = [dict(label=atlas_labels[k], voxels=int(counts[k]), percent=float(100*counts[k]/seed.sum()))
                     for k in np.argsort(counts)[::-1] if counts[k]]
    result = dict(sources=[identity(root/'masks'/name) for name in names], background=identity(args.background),
        seed_voxels=int(seed.sum()), seed_volume_ml=float(seed.sum()*voxel_volume/1000),
        centroid_mm=xyz.mean(0).tolist(), bounds_mm=[xyz.min(0).tolist(), xyz.max(0).tolist()],
        intersection_voxels=int(intersection.sum()), intersection_26_connected_clusters=count,
        intersection_matches_saved_task_masks=True, seed_is_exact_intersection_cluster=True,
        atlas=identity(atlas_path), atlas_label_file=identity(label_path), atlas_overlap=atlas_overlap,
        note='This verifies the tracked clusterthresh seed; historical wrappers also name a separate seed-AIns.nii.gz absent from this checkout.')
    (record/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    background = nib.as_closest_canonical(nib.load(args.background)); bg = background.get_fdata()
    atlas_display = resample_from_to(atlas, background, order=0).get_fdata() == insula_index
    shown = []
    for data in (seed, intersection & ~seed):
        image = nib.Nifti1Image(data.astype('uint8'), seed_image.affine)
        shown.append(resample_from_to(image, background, order=0).get_fdata() > .5)
    coords = tuple(np.round(xyz.mean(0)).astype(int)); ijk = np.round(nib.affines.apply_affine(np.linalg.inv(background.affine), coords)).astype(int)
    positions = [background.affine[k, 3]+np.arange(background.shape[k])*background.affine[k, k] for k in range(3)]
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2))
    vmax = np.percentile(bg[bg > 0], 99.5)
    for ax, axis, label in zip(axes, range(3), ('Sagittal', 'Coronal', 'Axial')):
        horizontal, vertical = ((1, 2), (0, 2), (0, 1))[axis]
        extent = (positions[horizontal][0]-.5, positions[horizontal][-1]+.5,
                  positions[vertical][0]-.5, positions[vertical][-1]+.5)
        plane = lambda a:np.take(a, ijk[axis], axis=axis).T
        ax.imshow(plane(bg), cmap='gray', vmin=0, vmax=vmax, origin='lower', extent=extent)
        for data, color in ((shown[1], '#21cbd7'), (shown[0], '#ffad20')):
            overlay = np.ma.masked_where(~plane(data), np.ones(plane(data).shape))
            ax.imshow(overlay, cmap=ListedColormap([color]), origin='lower', extent=extent, interpolation='nearest', alpha=.9, vmin=0, vmax=1)
        if plane(atlas_display).any():
            ax.contour(positions[horizontal], positions[vertical], plane(atlas_display), levels=[.5], colors=['#66f28a'], linewidths=1.15)
        ax.set_title(f'{label}: {"xyz"[axis]} = {coords[axis]} mm', fontsize=12)
        ax.set_xlim(((-100, 85), (-85, 85), (-85, 85))[axis]); ax.set_ylim(((-55, 90), (-55, 90), (-110, 85))[axis])
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
        if axis in (1, 2):
            ax.text(.03, .94, 'L', color='white', transform=ax.transAxes); ax.text(.94, .94, 'R', color='white', transform=ax.transAxes)
        else:
            ax.text(.03, .94, 'P', color='white', transform=ax.transAxes); ax.text(.94, .94, 'A', color='white', transform=ax.transAxes)
    fig.suptitle('Original task-intersection anterior-insula seed', fontsize=16, y=.97)
    fig.legend(handles=[Patch(color='#ffad20', label='Left insula seed: 127 voxels (2.75 mL)'),
                        Patch(color='#21cbd7', label='Other Trust ∩ UGR voxels'),
                        Line2D([0], [0], color='#66f28a', label='Anatomical insula (Harvard–Oxford 25%)')], loc='lower center', bbox_to_anchor=(.5, .075), ncol=3, frameon=False, fontsize=9)
    fig.text(.5, .025, 'Saved masks shown without smoothing; nearest-neighbor resampling for display on MNI152 T1. Neurological view: left is left.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .14, 1, .92)); fig.savefig(record/'original_insula_intersection.png', dpi=180); plt.close(fig)
    print(json.dumps({k:v for k,v in result.items() if k not in ('sources', 'background')}, indent=2))


if __name__ == '__main__': main()
