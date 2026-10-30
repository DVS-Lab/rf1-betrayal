#!/usr/bin/env python3
"""Show existing spatial maps and model-order correspondences without new ICA."""
from pathlib import Path
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np

import run_network_timecourse_analysis as n
from run_bold_controls import save


def main():
    source = n.RECORDS/'linux2-completed-ICA-exploration-20261029T185120Z-0fea45'
    output = n.RECORDS/'2026-10-29-component-interpretation'; output.mkdir(exist_ok=True)
    path = source/'primary_maps/joint_bold_bilateral_dim-5.nii.gz'
    image = nib.load(path); maps = image.get_fdata(); mask = np.any(maps != 0, axis=3)
    fsldir = Path(os.environ.get('FSLDIR', '/Users/tug87422/fsl'))
    background = fsldir/'data/standard/MNI152_T1_2mm_brain.nii.gz'
    bg = resample_from_to(nib.load(background), (maps.shape[:3], image.affine), order=1).get_fdata()
    coords = [image.affine[j, 3]+np.arange(maps.shape[j])*image.affine[j, j] for j in range(3)]
    vmax_bg = np.percentile(bg[bg > 0], 99)
    names = ('Anterior', 'Mid/anterior', 'Posterior/dorsal', 'Ventral mid/posterior', 'Central/inferior')
    fig, axes = plt.subplots(2, 5, figsize=(15, 7.2))
    for k, name in enumerate(names):
        v = maps[..., k]; normalized = v/np.max(abs(v))
        left = mask & (coords[0][:, None, None] < 0)
        ix = np.unravel_index(np.argmax(np.where(left, abs(v), -1)), v.shape)[0]
        peak = np.unravel_index(np.argmax(abs(v)), v.shape); iy = peak[1]
        for row, (data, anatomy, valid, horizontal, vertical) in enumerate((
            (normalized[ix, :, :].T, bg[ix, :, :].T, mask[ix, :, :].T, coords[1], coords[2]),
            (normalized[:, iy, :].T, bg[:, iy, :].T, mask[:, iy, :].T, coords[0], coords[2]))):
            ax = axes[row, k]; extent = (horizontal[0], horizontal[-1], vertical[0], vertical[-1])
            ax.imshow(anatomy, origin='lower', extent=extent, cmap='gray', vmin=0, vmax=vmax_bg)
            ax.contour(horizontal, vertical, valid.astype(float), [.5], colors='#73b8a5', linewidths=.6)
            shown = np.ma.masked_where(~valid | (abs(data) < .3), data)
            artist = ax.imshow(shown, origin='lower', extent=extent, cmap='RdBu_r', vmin=-1, vmax=1, alpha=.92)
            ax.set_xlim((-30, 40) if row == 0 else (-55, 55)); ax.set_ylim(-25, 30)
            ax.set(xticks=[], yticks=[])
            if row == 0: ax.set_title(f'IC{k+1}\n{name}', fontsize=11)
            ax.set_xlabel(f'x={coords[0][ix]:.0f} mm; posterior to anterior' if row == 0 else f'y={coords[1][iy]:.1f} mm; left to right', fontsize=8)
    axes[0, 0].set_ylabel('Sagittal slice at left loading peak\nSuperior up', fontsize=10)
    axes[1, 0].set_ylabel('Coronal slice at loading peak\nSuperior up', fontsize=10)
    fig.suptitle('Bilateral five-component insula ICA: spatial patterns, not five disjoint parcels', fontsize=15)
    fig.text(.5, .035, 'Colors: signed loading / largest absolute loading within that component. Display |loading| >= 30% of peak; descriptive, not a significance threshold.', ha='center', fontsize=8)
    fig.text(.5, .012, 'Green contour: defining insula support. Standard MNI T1 underlay. Positive/negative lobes can overlap across components.', ha='center', fontsize=8)
    fig.subplots_adjust(left=.065, right=.93, bottom=.12, top=.84, hspace=.32, wspace=.12)
    cax = fig.add_axes([.95, .23, .012, .48]); fig.colorbar(artist, cax=cax, ticks=[-1, 0, 1], label='Relative signed loading')
    fig.savefig(output/'bilateral_five_component_anatomy.png', dpi=170); plt.close(fig)

    comparisons = json.loads((source/'fixed_dimension_comparisons.json').read_text())
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios':[3, 7]})
    for ax, dim in zip(axes, (3, 7)):
        item = next(r for r in comparisons if r['mode'] == 'joint_bold' and r['roi'] == 'bilateral' and r['sensitivity_dimension'] == dim)
        corr = np.array(item['all_pair_spatial_r']); artist = ax.imshow(abs(corr), vmin=0, vmax=1, cmap='Blues', aspect='equal')
        for i in range(5):
            for j in range(dim): ax.text(j, i, f'{corr[i,j]:+.2f}', ha='center', va='center', fontsize=9,
                                        color='white' if abs(corr[i,j])>.65 else 'black')
        for match in item['matches']:
            ax.add_patch(Rectangle((match['sensitivity_component']-1.5, match['primary_component']-1.5), 1, 1,
                                   fill=False, edgecolor='#d57222', lw=2))
        ax.set(xticks=range(dim), xticklabels=['IC'+str(k) for k in range(1, dim+1)], yticks=range(5),
               yticklabels=['5D IC'+str(k) for k in range(1, 6)], xlabel=f'{dim}-component solution', title=f'Five versus {dim} components')
    fig.suptitle('Same anatomy, changing model order: match spatial patterns rather than IC numbers', fontsize=13)
    fig.text(.5, .025, 'Cell text: signed spatial r. Color: |r|. Orange outlines: one-to-one matches. Unmatched lower-order components can be mixtures.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .07, .95, .92)); cax = fig.add_axes([.96, .22, .012, .53]); fig.colorbar(artist, cax=cax, label='Absolute spatial r')
    fig.savefig(output/'bilateral_dimension_correspondence.png', dpi=170); plt.close(fig)
    save(output/'figure-provenance.json', dict(map=n.g.core.identity(path), comparison=n.g.core.identity(source/'fixed_dimension_comparisons.json'),
        background=n.g.core.identity(background), map_threshold_fraction=.3, physical_map_data_changed=False,
        task_tests_at_other_dimensions='Not yet performed with the current pooled standardized OLS model.'))
    print(output)


if __name__ == '__main__': main()
