#!/usr/bin/env python3
"""Collect original manuscript group maps, live seeds and UGR design diagnostics."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
import traceback
import uuid

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'code/insula_ica'))
import insula_ica as core
from run_bold_controls import save, publish


def contrast_diagnostics(design, contrasts):
    """Audit original raw-unit contrasts after scaling columns for numerical SVD."""
    core.require(np.isfinite(design).all() and np.isfinite(contrasts).all(), 'Nonfinite model or contrast')
    centered = design-design.mean(axis=0)
    sd = centered.std(axis=0, ddof=1)
    active = sd > 1e-12*np.maximum(1., np.max(np.abs(design), axis=0))
    x = centered[:, active]/sd[active]
    core.require(x.shape[1] > 0, 'No nonconstant model columns')
    # Original beta = standardized-design beta / SD(X); do NOT silently use
    # unit standardized-contrast weights for the manuscript's original contrast.
    cs = contrasts[:, active]/sd[active]
    u, singular, vt = np.linalg.svd(x, full_matrices=False)
    results = []
    for cutoff in (1e-10, 1e-8, 1e-6):
        kept = singular > singular[0]*cutoff; basis = vt[kept]
        projected = (cs@basis.T)@basis
        residual = np.linalg.norm(cs-projected, axis=1)/np.maximum(np.linalg.norm(cs, axis=1), 1e-15)
        inactive_weights = np.any(np.abs(contrasts[:, ~active]) > 1e-10, axis=1)
        for k, value in enumerate(residual):
            results.append(dict(contrast=k+1, cutoff=cutoff, rank=int(kept.sum()), active_columns=int(active.sum()),
                relative_nonestimable_fraction=float(value), weights_on_constant_columns=bool(inactive_weights[k]),
                requires_review=bool(value > 1e-5 or inactive_weights[k])))
    return results


def collect(record):
    manifest = []; maps = record/'group_maps'; sources = record/'sources'
    maps.mkdir(); sources.mkdir()

    def capture(path, role):
        path = Path(path)
        entry = dict(source=str(path), role=role, exists=path.is_file())
        if path.is_file():
            try:
                if path.name.endswith(('.nii', '.nii.gz')):
                    im = nib.load(path)
                    core.require(len(im.shape) == 3, 'Only group/seed 3D imaging may be copied')
                    entry.update(shape=list(im.shape), affine=im.affine.tolist())
                    directory = maps
                else: directory = sources
                identity = core.identity(path); tag = identity['sha256'][:12]
                target = directory/(tag+'-'+path.name)
                if not target.exists(): shutil.copyfile(path, target)
                entry.update(identity=identity, copy=str(target.relative_to(record)))
            except Exception as exc: entry['error'] = str(exc)
        manifest.append(entry)
        return entry

    for folder in (ROOT/'masks', ROOT/'masks/archive'):
        if folder.is_dir():
            for path in sorted(folder.iterdir()):
                if path.is_file() and path.name.endswith(('.nii', '.nii.gz')) and any(term in path.name.lower() for term in ('insula', 'ains', 'intersection', 'overlap_figure')):
                    capture(path, 'live_seed_or_overlap')

    group_root = ROOT/'derivatives/fsl/L3-act'; cases = []; group_paths = set()
    # Bounded to the actual manuscript group root; no scan of other projects.
    for task, model, copes in (('ugr', '3', (11, 12)), ('trust', '01', (10,))):
        for parent in sorted(group_root.glob(f'L3_model-{model}_task-{task}_type-act-n132-cov-*')):
            if not parent.is_dir(): continue
            for group in sorted(parent.glob('*.gfeat')):
                if not any(f'-cope-{cope}_cname-' in group.name for cope in copes): continue
                group_paths.add(group); case = dict(path=str(group), task=task)
                for folder in (group, group/'cope1.feat'):
                    for name in ('design.fsf', 'design.mat', 'design.con', 'design.grp'):
                        if (folder/name).is_file(): capture(folder/name, 'group_design')
                fsf = next((f/'design.fsf' for f in (group, group/'cope1.feat') if (f/'design.fsf').is_file()), group/'design.fsf')
                if fsf.is_file():
                    settings = core.fsl_settings(fsf)
                    case['saved_input_paths'] = {k:v for k,v in settings.items() if k.startswith('feat_files(')}
                    case['saved_contrast_names'] = {k:v for k,v in settings.items() if k.startswith('fmri(conname_')}
                    case['thresholds'] = {k:v for k,v in settings.items() if k in ('fmri(thresh)', 'fmri(prob_thresh)', 'fmri(z_thresh)', 'fmri(mixed_yn)')}
                for name in ('mask.nii.gz', 'cluster_mask_zstat1.nii.gz', 'cluster_mask_zstat2.nii.gz',
                             'thresh_zstat1.nii.gz', 'thresh_zstat2.nii.gz', 'stats/zstat1.nii.gz', 'stats/zstat2.nii.gz'):
                    capture(group/'cope1.feat'/name, 'original_group_map')
                cases.append(case)
    # Record the exact script-named intersection parents even when missing.
    script = ROOT/'code/get_ugr-trust_intersection.sh'; capture(script, 'intersection_script')
    import re
    expected = []
    for value in re.findall(r'^(?:ug|trust)=\$maindir(/[^\n]+)', script.read_text(), re.MULTILINE):
        path = ROOT/value.lstrip('/'); expected.append(capture(path, 'script_named_intersection_parent'))
    save(record/'group_cases.json', cases)

    registry_path = ROOT/'code/insula_ica/records/linux2-network-timecourses-20261030T012523Z-613cec/accepted-inputs.json'
    registry = json.loads(registry_path.read_text()); capture(registry_path, 'manuscript_run_registry')
    rows = [r for r in registry if r['task'] == 'ugr']
    core.require(len(rows) == 249 and len({r['subject'] for r in rows}) == 132, 'Original accepted UGR cohort changed')
    audits = []; failures = []; timing = []
    for row in rows:
        subject = str(row['subject']).removeprefix('sub-')
        feat = ROOT/f'derivatives/fsl/sub-{subject}/ses-01/L1_task-ugr_model-3_type-act_run-{row["run"]}_sm-5.feat'
        try:
            settings = core.fsl_settings(feat/'design.fsf'); design = core.fsl_matrix(feat/'design.mat'); contrasts = core.fsl_matrix(feat/'design.con')
            core.require(design.shape[1] == contrasts.shape[1], 'Design/contrast dimensions differ')
            core.require(design.shape[0] == row['nvol'], 'Original model volume count differs from accepted run')
            timing.append(dict(key=row['key'], expected_acquisition_TR=row['tr'], historical_model_TR=float(settings['fmri(tr)']),
                               mismatch=not np.isclose(float(settings['fmri(tr)']), float(row['tr']), atol=1e-5)))
            for item in contrast_diagnostics(design, contrasts):
                audits.append(dict(key=row['key'], contrast_name=settings.get(f'fmri(conname_real.{item["contrast"]})', ''), **item))
            if subject in ('10817', '11031'):
                for name in ('design.fsf', 'design.mat', 'design.con', 'confoundevs.txt'):
                    capture(feat/name, 'affected_original_L1_model')
                for k, value in settings.items():
                    if re.fullmatch(r'fmri\(custom\d+\)', k): capture(value, 'affected_original_event_file')
        except Exception as exc: failures.append(dict(key=row['key'], path=str(feat), error=str(exc)))
    core.write_table(record/'original_UGR_contrast_audit.tsv', audits)
    save(record/'original_UGR_timing.json', timing); save(record/'missing_or_failed_L1.json', failures)
    save(record/'manifest.json', manifest)
    flagged = sorted({r['key'] for r in audits if r['requires_review']})
    state = dict(state='original_group_input_audit_complete_review_required', UGR_runs_expected=249,
        UGR_runs_audited=len({r['key'] for r in audits}), UGR_runs_failed=len(failures), flagged_design_runs=flagged,
        TR_mismatches=[r for r in timing if r['mismatch']], group_cases=len(cases),
        script_named_parents_found=sum(r['exists'] for r in expected), script_named_parents_expected=len(expected),
        copied_files=sum('copy' in r for r in manifest), collection_errors=[r for r in manifest if r.get('error')],
        data_or_analysis_modified=False, FEAT_or_group_jobs_launched=False)
    save(record/'state.json', state)
    (record/'README.md').write_text('# Original manuscript group-map and UGR audit\n\n'
        'Read-only inventory of live seeds, original UGR/Trust N132 activation group maps (UGR COPE11/12, Trust COPE10), '
        'saved group designs, and all 249 original sm-5 UGR designs. See manifest.json for original paths, hashes and copies. '
        'Map names are not assumed to establish effect direction; group_cases.json preserves saved contrast names and input paths. '
        'Only 3D group/seed images and small model files are copied. No subject 4D images, FEAT jobs, or existing output changes.\n\n'
        'The original contrast estimability diagnostic scales design columns for numerical SVD and transforms original contrast weights by the same SDs. '
        'It reports three numerical cutoffs, including a conservative 1e-6 diagnostic, rather than changing original models or substituting standardized task contrasts. '
        'Constant-column weights and near-dependent contrasts require review. A review flag is not a new exclusion. '
        'Timing compares the original model with the accepted corrected acquisition-TR registry. '
        'Missing originals are recorded and never replaced silently by the newer sm-0 ICA stream.\n')
    print(json.dumps(state, indent=2), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--publish', action='store_true'); args = p.parse_args()
    records = ROOT/'code/ugr_reanalysis/records'; records.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    record = records/('linux2-original-group-audit-'+stamp); record.mkdir()
    try: collect(record)
    except Exception as exc:
        save(record/'state.json', dict(state='original_group_audit_failed', error=str(exc)))
        (record/'failure.txt').write_text(traceback.format_exc()); print('Audit failed:', exc, flush=True)
    if args.publish:
        try: publish(record, 'Collect original seed/group maps and audit manuscript UGR design contrasts')
        except Exception as exc: (record/'publication-error.txt').write_text(str(exc)); print('Git publication failed:', exc, flush=True)
    print('Record:', record, flush=True)
    raise SystemExit(0 if json.loads((record/'state.json').read_text())['state'].endswith('review_required') else 1)


if __name__ == '__main__': main()
