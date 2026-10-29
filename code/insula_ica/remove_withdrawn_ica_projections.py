#!/usr/bin/env python3
"""Remove only the withdrawn ICA task-projection artifacts on linux2.

Git pull removes the tracked tables/figures. This also removes the same specific
outputs in ignored exploration folders and strips obsolete estimates from run
completion JSON. Never removes any FEAT COPE, ICA map, residual, or timecourse.
"""
import json
from pathlib import Path
import subprocess

import insula_ica as core
from run_bold_controls import save

NAMES = ('run_task_projections.tsv', 'participant_task_projections.tsv',
         'group_task_projections.tsv', 'task_cope_provenance.json',
         'figures/joint_bold_task_projections.png', 'figures/native_film_task_projections.png')


def clean():
    processes = subprocess.check_output(['ps', '-eo', 'args'], text=True)
    core.require(not any('explore_completed_ica.py' in line and 'python' in line for line in processes.splitlines()),
                 'An exploration process is active; wait for it to finish before cleanup')
    removed = []; rewritten = []
    records = core.ROOT/'code/insula_ica/records'
    for record in records.glob('linux2-completed-ICA-exploration-*'):
        for name in NAMES:
            path = record/name
            if path.is_file(): path.unlink(); removed.append(str(path))
        state_path = record/'state.json'
        if not state_path.is_file(): continue
        state = json.loads(state_path.read_text()); output = Path(state['output_root']).resolve()
        core.require(output.is_relative_to((core.ROOT/'derivatives/insula_ica').resolve()) and
                     output.name.startswith('completed-ICA-exploration-'), 'Unexpected exploration root; no cleanup')
        for name in NAMES:
            path = output/name
            if path.is_file(): path.unlink(); removed.append(str(path))
        for marker in (output/'runs').glob('*/complete.json'):
            value = json.loads(marker.read_text())
            changed = value.get('result', {}).pop('effects', None) is not None
            changed = (value.get('signature', {}).pop('cope', None) is not None) or changed
            if changed: save(marker, value); rewritten.append(str(marker))
    print(f'Withdrawn task outputs removed: {len(removed)}; run markers cleaned: {len(rewritten)}')
    return dict(removed=removed, completion_markers_cleaned=rewritten,
                retained='All FEAT maps, original data, ICA fits, connectivity maps, residuals and timecourses')


if __name__ == '__main__':
    result = clean()
    output = core.ROOT/'code/insula_ica/output/withdrawn-task-cleanup.json'
    output.parent.mkdir(exist_ok=True); save(output, result)
    print('Cleanup record:', output)
