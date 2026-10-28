#!/usr/bin/env python3
"""Export review images for flagged Trust runs, without fitting or excluding runs."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import traceback
import uuid

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

import insula_ica as core
from run_bold_controls import publish, save


SOURCE = 'linux2-trust-regeneration-20261027T162428Z-ef53c5'


def montage(qc, destination, fsldir):
    reference = nib.load(qc['historical_frame']['reference'])
    bold = nib.load(qc['historical_frame']['bold'])
    old_mask = nib.load(Path(qc['historical_frame']['reference']).parent / 'mask.nii.gz')
    new_mask = nib.load(Path(qc['residual']).parents[1] / 'mask.nii.gz')
    core.require(all(core.same_grid(reference, im) for im in (bold, old_mask, new_mask)),
                 'Review requires matching image grids: ' + qc['key'])
    # Preserve the original orientation for overlays; reject unexpected obliquity.
    affine = reference.affine
    core.require(np.allclose(affine[:3, :3], np.diag(np.diag(affine[:3, :3])), atol=1e-5),
                 'Unexpected oblique grid; inspect in an image viewer: ' + qc['key'])
    core.require(all(np.diag(affine)[:3] > 0), 'Expected increasing RAS coordinates')
    rois, _, _ = core.anatomical_masks(
        fsldir / 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz',
        fsldir / 'data/atlases/HarvardOxford-Cortical.xml', reference)
    insula = rois['bilateral']
    masks = [old_mask.get_fdata() > 0, new_mask.get_fdata() > 0]
    index = qc['historical_frame']['frame_index']
    images = [reference.get_fdata(), np.asarray(bold.dataobj[..., index], dtype=float)]
    z_mm = [-15, -5, 5, 15, 25, 35]
    slices = [int(round((z-affine[2, 3])/affine[2, 2])) for z in z_mm]
    core.require(all(0 <= z < reference.shape[2] for z in slices), 'Requested slices outside field of view')
    fig, axes = plt.subplots(2, 6, figsize=(15, 6), facecolor='black')
    for row, (data, mask) in enumerate(zip(images, masks)):
        inside = data[masks[0] & np.isfinite(data)]
        core.require(len(inside) > 0, 'No finite review-image intensities')
        low, high = np.percentile(inside, [1, 99])
        for col, (z, mm) in enumerate(zip(slices, z_mm)):
            ax = axes[row, col]
            ax.imshow(data[:, :, z].T, origin='lower', cmap='gray', vmin=low, vmax=max(high, low+1))
            for overlay, color in ((mask, '#ffae00'), (insula, '#00ffff')):
                plane = overlay[:, :, z].T
                if plane.any() and not plane.all():
                    ax.contour(plane, levels=[.5], colors=[color], linewidths=.55)
            ax.axis('off')
            if row == 0:
                ax.set_title(f'z = {mm} mm', color='white', fontsize=10)
            if col == 0:
                ax.text(.02, .02, 'Historical reference' if row == 0 else f'Regenerated frame {index}',
                        transform=ax.transAxes, color='white', fontsize=9,
                        bbox=dict(facecolor='black', alpha=.7, edgecolor='none'))
    cov = qc['coverage']
    fig.suptitle(qc['key'] + '\n' +
        f"Historical spatial r={qc['historical_frame']['spatial_r']:.3f}; insula coverage " +
        ', '.join(f'{name} {100*cov[name]["retained_fraction"]:.1f}%' for name in ('left', 'right', 'bilateral')),
        color='white', fontsize=12)
    fig.text(.5, .025, 'Orange: corresponding FEAT mask | Cyan: anatomical insula | Left is image left (neurological)\n'
             'Each row has its own intensity window. Historical differences alone do not establish a registration failure.',
             ha='center', color='white', fontsize=9)
    fig.subplots_adjust(left=.01, right=.99, top=.82, bottom=.10, wspace=.02, hspace=.04)
    fig.savefig(destination, dpi=110, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    source = core.ROOT / 'code/insula_ica/records' / SOURCE
    qcs = [json.loads(p.read_text()) for p in sorted((source / 'qc').glob('*.json'))]
    core.require(len(qcs) == 241, 'Expected all 241 source QC reports')
    selected = [q for q in qcs if q['review_flags']]
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    record = core.ROOT / 'code/insula_ica/records' / ('linux2-trust-visual-qc-' + stamp)
    record.mkdir(parents=True)
    fsldir = Path(os.environ.get('FSLDIR', '/usr/local/fsl'))
    results, errors = [], []
    for qc in sorted(selected, key=lambda q: (q['coverage']['bilateral']['retained_fraction'], q['key'])):
        item = dict(key=qc['key'], source_qc=str(source / 'qc' / (qc['key'] + '.json')),
                    review_flags=qc['review_flags'], fmriprep_report=qc['fmriprep_report'])
        try:
            image = qc['key'] + '.png'
            montage(qc, record / image, fsldir)
            item.update(image=image, state='exported')
        except Exception as exc:
            item.update(state='failed', error=str(exc))
            error_file = record / (qc['key'] + '.error.txt')
            error_file.write_text(traceback.format_exc())
            errors.append(item)
        results.append(item)
    save(record / 'index.json', dict(source_record=str(source), runs=results,
                                    review_state='not_yet_visually_reviewed', exclusions=[]))
    lines = ['# Trust images for visual QC', '',
             'Review bundle only. No preprocessing, FEAT, ICA, or exclusions were performed.', '',
             'The two sub-10834 runs have low insula coverage; inspect these first. '
             'Other flags compare against the historical reference and need visual interpretation.', '',
             'The images compare the retained historical reference with the corresponding regenerated BOLD frame. '
             'Orange contours are each row\'s FEAT mask; cyan contours are the anatomical insula. '
             'Intensity windows are independent. No visual pass/fail decisions have been assigned.', '']
    for item in results:
        lines += ['## ' + item['key'], '', ', '.join(item['review_flags']), '',
                  f"![Review montage]({item['image']})" if item['state'] == 'exported' else 'EXPORT ERROR: ' + item['error'], '',
                  'Full fMRIPrep report on linux2: `' + item['fmriprep_report'] + '`', '']
    (record / 'README.md').write_text('\n'.join(lines))
    print(f'Exported {len(selected)-len(errors)}/{len(selected)} review images; errors: {len(errors)}', flush=True)
    print('Record:', record, flush=True)
    if args.publish:
        publish(record, 'Export flagged Trust runs for visual QC')
    raise SystemExit(bool(errors))


if __name__ == '__main__':
    main()
