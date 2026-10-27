"""Numerical comparisons for independent, same-grid BOLD controls."""
from pathlib import Path

import nibabel as nib
import numpy as np

from diagnose_film_replay import agreement
from rebuild_trust import require


def same_grid(a, b):
    return a.shape[:3] == b.shape[:3] and np.allclose(a.affine, b.affine, atol=1e-5)


def geometry(image):
    return dict(shape=list(image.shape), zooms=[float(x) for x in image.header.get_zooms()],
                affine=image.affine.tolist(), stored_dtype=str(image.get_data_dtype()))


def frame_metrics(reference, bold, mask, index):
    a, b, m = nib.load(reference), nib.load(bold), nib.load(mask)
    info = dict(reference=str(reference), bold=str(bold), frame_index=index,
                same_grid=bool(same_grid(a, b) and same_grid(a, m)),
                reference_geometry=geometry(a), bold_geometry=geometry(b))
    if not info['same_grid']:
        return info
    require(len(b.shape) == 4 and index < b.shape[3], 'Incompatible BOLD dimensions')
    x, y = a.get_fdata().squeeze(), np.asarray(b.dataobj[..., index], dtype=float)
    valid = (m.get_fdata() > 0) & np.isfinite(x) & np.isfinite(y)
    info.update(agreement(x[valid], y[valid]))
    den = float(x[valid] @ x[valid])
    if den:
        scale = float(x[valid] @ y[valid] / den)
        info.update(global_scale=scale, relative_rms_after_global_scale=float(
            np.linalg.norm(y[valid]-scale*x[valid]) / np.sqrt(den)))
    return info


def series_metrics(left, right, mask):
    """Compare all time points, with voxel-demeaned and temporal-correlation checks.

    Read each compressed series once, then accumulate in small spatial chunks.
    Do not resample, normalize, filter, or fit another nuisance model.
    """
    a, b, m = nib.load(left), nib.load(right), nib.load(mask)
    info = dict(left=str(left), right=str(right), mask=str(mask),
                left_shape=list(a.shape), right_shape=list(b.shape),
                same_grid=bool(same_grid(a, b) and same_grid(a, m)),
                left_geometry=geometry(a), right_geometry=geometry(b))
    if not info['same_grid'] or a.shape != b.shape or len(a.shape) != 4:
        info['compared'] = False
        return info
    # Preserve native precision; float32 conversion could hide small differences.
    aa, bb = np.asanyarray(a.dataobj), np.asanyarray(b.dataobj)
    positions = np.flatnonzero(m.get_fdata().reshape(-1) > 0)
    aa, bb = aa.reshape(-1, a.shape[-1]), bb.reshape(-1, b.shape[-1])
    totals = dict(xx=0., yy=0., xy=0., dd=0., centered_xx=0., centered_dd=0.)
    correlations, max_abs, equal, n, skipped = [], 0., True, 0, 0
    for start in range(0, len(positions), 2048):
        idx = positions[start:start+2048]
        x, y = aa[idx].astype(float), bb[idx].astype(float)
        finite = np.isfinite(x).all(axis=1) & np.isfinite(y).all(axis=1)
        skipped += int((~finite).sum())
        x, y = x[finite], y[finite]
        if not len(x):
            continue
        n += len(x)
        diff = y-x
        equal = equal and bool(np.array_equal(x, y))
        max_abs = max(max_abs, float(np.max(np.abs(diff))))
        for name, values in (('xx', x*x), ('yy', y*y), ('xy', x*y), ('dd', diff*diff)):
            totals[name] += float(values.sum())
        x -= x.mean(axis=1, keepdims=True)
        y -= y.mean(axis=1, keepdims=True)
        totals['centered_xx'] += float((x*x).sum())
        totals['centered_dd'] += float(((y-x)**2).sum())
        den = np.sqrt((x*x).sum(axis=1) * (y*y).sum(axis=1))
        valid = den > 0
        correlations.extend(((x*y).sum(axis=1)[valid] / den[valid]).tolist())
    require(n > 0, 'No finite voxels in series comparison')
    scale = totals['xy']/totals['xx'] if totals['xx'] else None
    scaled = max(0., totals['yy']-2*scale*totals['xy']+scale*scale*totals['xx']) if scale is not None else None
    info.update(compared=True, voxels=n, nonfinite_voxels_skipped=skipped,
                exact_numeric_equality=equal, max_absolute_difference=max_abs,
                relative_rms_difference=np.sqrt(totals['dd']/totals['xx']) if totals['xx'] else None,
                global_scale=scale,
                relative_rms_after_global_scale=np.sqrt(scaled/totals['xx']) if scaled is not None else None,
                voxel_demeaned_relative_rms=np.sqrt(totals['centered_dd']/totals['centered_xx']) if totals['centered_xx'] else None,
                voxel_temporal_correlation_percentiles=dict(zip(('p05', 'median', 'p95'),
                    np.percentile(correlations, [5, 50, 95]).tolist())) if correlations else None)
    return info


def find_bold(output, row):
    expected = Path(output) / 'fmriprep' / row['subject'] / 'ses-01/func' / Path(row['new_bold']).name
    if expected.is_file():
        return expected
    options = [p for p in (Path(output) / 'fmriprep' / row['subject']).rglob(
        row['key'] + '_*space-MNI152NLin6Asym*_desc-preproc_bold.nii.gz') if '_echo-' not in p.name]
    require(len(options) == 1, f'Expected one combined MNI BOLD for {row["key"]}; found {len(options)}')
    return options[0]
