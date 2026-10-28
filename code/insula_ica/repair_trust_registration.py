#!/usr/bin/env python3
"""Test sub-10834 with a BOLD-derived reference instead of the SBRef.

One participant, both accepted runs, fresh output/work directories. This is a
candidate correction pending visual review, not an automatic replacement.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import traceback
import uuid

import insula_ica as core
import rebuild_trust as rebuild
from diagnose_trust_registration import collect, sbref_selection, SOURCE, SUBJECT
from export_trust_qc import montage
from regenerate_trust import production_command, qc_run
from run_bold_controls import publish, save


WORK_BASE = Path('/ZPOOL/data/scratch') / os.environ.get('USER', 'tug87422') / 'rf1-betrayal-registration'


def candidate_command(runtime, image, upstream, output, work):
    return production_command(runtime, image, upstream, output, work, SUBJECT) + ['--ignore', 'sbref']


def candidate_rows(rows, output):
    selected = [dict(r) for r in rows if r['subject'] == SUBJECT]
    core.require(len(selected) == 2 and {int(r['run']) for r in selected} == {1, 2},
                 'Expected both accepted sub-10834 runs')
    for row in selected:
        # Preserve the original BIDS basename, scientific model and nuisance path.
        row['new_bold'] = str(output / 'fmriprep-24' / SUBJECT / 'ses-01/func' / Path(row['new_bold']).name)
    return selected


def run(record, state, stamp):
    source = core.ROOT / 'code/insula_ica/records' / SOURCE
    all_rows = json.loads((source / 'input-registry.json').read_text())
    original = [r for r in all_rows if r['subject'] == SUBJECT]
    save(record / 'cohort-sbref-selection-audit.json', sbref_selection(all_rows))
    baseline = record / 'baseline-registration-reports'; baseline.mkdir()
    baseline_assets, baseline_errors = collect(source, baseline)
    core.require(baseline_assets > 0 and not baseline_errors, 'Baseline registration report export incomplete')
    output = core.ROOT / 'derivatives/insula_ica' / ('registration-' + SUBJECT + '-ignore-sbref-' + stamp)
    work = WORK_BASE / output.name
    core.require(not output.exists() and not work.exists(), 'Correction requires fresh output/work directories')
    rows = candidate_rows(original, output)
    core.require(all(r['warpkit']['ready'] for r in rows), 'Expected verified existing Warpkit inputs')
    for row in rows:
        core.require(core.sha(row['confounds']) == row['confounds_sha256'], 'Original preserved confounds changed')
    before = {r['new_bold']: core.identity(r['new_bold'], content=False) for r in original}
    save(record / 'original-output-identities.json', before)
    output.mkdir(parents=True); work.mkdir(parents=True)
    (output / 'filters').mkdir()
    save(output / 'filters' / (SUBJECT + '.json'), rebuild.filter_for(rows))
    runtime = shutil.which('apptainer') or shutil.which('singularity')
    core.require(runtime and shutil.which('feat') and os.environ.get('FSLDIR'), 'Normal container and FSL setup required')
    image = Path('/ZPOOL/data/tools/fmriprep-24.1.1.simg')
    version = subprocess.check_output([runtime, 'run', '--cleanenv', str(image), '--version'],
                                     text=True, stderr=subprocess.STDOUT)
    rebuild.record_container(output, image, version)
    container = json.loads((output / 'container.json').read_text())
    expected = json.loads((source / 'container.json').read_text())
    core.require(container['sha256'] == expected['sha256'], 'Candidate container differs from completed cohort')
    save(record / 'container.json', container)
    argv = candidate_command(runtime, image, Path('/ZPOOL/data/projects/rf1-sra-linux2'), output, work)
    plan = dict(subject=SUBJECT, work=str(work), argv=argv, keys=[r['key'] for r in rows])
    save(record / 'plan.json', [plan]); save(record / 'input-registry.json', rows)
    save(record / 'hypothesis.json', dict(
        diagnosis='Gross BOLD-to-T1w mismatch, poor BOLD-reference masking in both runs.',
        candidate='Ignore SBRef and estimate the reference from the BOLD series.',
        additional_filter_correction='SBRef queries now inherit each BOLD run instead of overriding run with a list.',
        changed_scientific_setting=['--ignore sbref'], all_other_command_settings_from='production_command',
        source_record=str(source), original_confounds_required=True,
        automatic_adoption=False, participant_exclusions=[],
        documentation='https://fmriprep.org/en/24.1.1/usage.html#workflow-configuration'))
    subprocess.run([runtime, 'exec', '--cleanenv', '-B', f'{output}:/rebuild:ro', str(image),
                    'python', '-c', rebuild.FILTER_PREFLIGHT, f'/rebuild/filters/{SUBJECT}.json'],
                   stdout=subprocess.DEVNULL, check=True)
    env = dict(os.environ)
    for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
        env.update({prefix+'TEMPLATEFLOW_HOME': '/opt/templateflow', prefix+'MPLCONFIGDIR': '/opt/mplconfigdir'})
    state.update(state='running_fmriprep', output_root=str(output), work=str(work))
    save(record / 'state.json', state)
    log = output / (SUBJECT + '.fmriprep.log')
    print('Running candidate correction:', SUBJECT, '| log:', log, flush=True)
    with log.open('w') as stream:
        subprocess.run(argv, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
    rebuild.verify_outputs(rows)
    diagnostic = record / 'registration-reports'; diagnostic.mkdir()
    assets, header_errors = collect(record, diagnostic)
    core.require(assets > 0 and not header_errors, 'Candidate registration report export incomplete')
    qc_dir = record / 'qc'; qc_dir.mkdir()
    comparisons = []
    for row in rows:
        state.update(state='fitting_original_joint_model', current_run=row['key'])
        save(record / 'state.json', state)
        rebuild.fit_run(row, output, strict_historical=False)
        qc = qc_run(row, output)
        save(qc_dir / (row['key'] + '.json'), qc)
        montage(qc, qc_dir / (row['key'] + '.png'), Path(os.environ['FSLDIR']))
        old = json.loads((source / 'qc' / (row['key'] + '.json')).read_text())
        comparisons.append(dict(key=row['key'], before_coverage=old['coverage'], candidate_coverage=qc['coverage'],
            before_spatial_r=old['historical_frame']['spatial_r'], candidate_spatial_r=qc['historical_frame']['spatial_r'],
            confounds_exact=qc['confounds']['exact_numeric_equality'],
            hemisphere_coverage_at_least_80_percent=all(qc['coverage'][r]['retained_fraction'] >= .8 for r in ('left','right'))))
    core.require(before == {r['new_bold']: core.identity(r['new_bold'], content=False) for r in original},
                 'Original completed BOLD identities changed during candidate run')
    save(record / 'comparison.json', comparisons)
    coverage_ok = all(x['hemisphere_coverage_at_least_80_percent'] for x in comparisons)
    state.update(state='candidate_ready_for_visual_review' if coverage_ok else 'candidate_still_fails_coverage',
                 finished_utc=datetime.now(timezone.utc).isoformat(), original_outputs_unchanged=True,
                 automatic_adoption=False)
    lines = ['# sub-10834 candidate correction: ignore SBRef', '',
        'Same fMRIPrep 24.1.1 image and production settings, with `--ignore sbref`. '
        'Both Trust runs were refitted with their original joint task+nuisance models. '
        'Candidate outputs remain separate; none have replaced the completed cohort.', '',
        '| Run | Before bilateral coverage | Candidate bilateral coverage | Before frame r | Candidate frame r |',
        '| --- | ---: | ---: | ---: | ---: |']
    for x in comparisons:
        lines.append(f"| {x['key']} | {100*x['before_coverage']['bilateral']['retained_fraction']:.1f}% | "
                     f"{100*x['candidate_coverage']['bilateral']['retained_fraction']:.1f}% | "
                     f"{x['before_spatial_r']:.3f} | {x['candidate_spatial_r']:.3f} |")
    lines += ['', 'Full left/right/bilateral values and exact-confound checks: comparison.json. '
              'Inspect both QC montages, baseline-registration-reports/ and registration-reports/ before adopting any candidate. '
              'Coverage or historical correlation alone does not establish successful registration.', '']
    for row in rows: lines += ['!['+row['key']+'](qc/'+row['key']+'.png)', '']
    (record / 'README.md').write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    record = core.ROOT / 'code/insula_ica/records' / ('linux2-registration-fix-' + SUBJECT + '-' + stamp)
    record.mkdir(parents=True)
    state = dict(state='preflight', started_utc=datetime.now(timezone.utc).isoformat(),
                 code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=core.ROOT,text=True).strip())
    save(record / 'state.json', state)
    try:
        run(record, state, stamp)
    except Exception as exc:
        state.update(state='failed', error=str(exc))
        (record / 'failure.txt').write_text(traceback.format_exc())
        print('ERROR:', exc, flush=True)
    save(record / 'state.json', state)
    if args.publish:
        try:
            publish(record, 'Test sub-10834 BOLD-derived reference correction')
        except Exception as exc:
            (record / 'publication-error.txt').write_text(str(exc))
            print('Publication failed; record retained:', exc, flush=True)
    print('Final state:', state['state'], '| record:', record, flush=True)
    raise SystemExit(0 if state['state'] == 'candidate_ready_for_visual_review' else 1)


if __name__ == '__main__':
    main()
