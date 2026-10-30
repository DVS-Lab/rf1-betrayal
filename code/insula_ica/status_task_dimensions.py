#!/usr/bin/env python3
"""Bounded 3/7-component task and insula-mean status without imaging reads."""
import argparse
from pathlib import Path
import status_insula_streams as s


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--brief', action='store_true'); p.parse_args()
    choices = sorted((s.ROOT/'code/insula_ica/records').glob('linux2-task-dimensions-*'))
    if not choices: print('No 3/7-component task sensitivity record found.'); return
    record = choices[-1]; state = s.read(record/'state.json'); rows = s.processes()
    runners, _, _ = s.pipeline(rows or [], ('task_dimension_sensitivity.py',), Path(state['output_root']))
    print('Record:', record.name)
    print('Runner:', 'unknown' if rows is None else ', '.join(f'PID {r["pid"]}, elapsed {r["elapsed"]}' for r in runners) if runners else 'not running')
    print('State:', state['state'], '| pilot passed:', state['pilot_passed'])
    print('Runs complete:', sum(r['state'] == 'complete' for r in state['runs'].values()), '/ 490')
    print('Runs reusing extraction:', sum(r.get('extraction_reused', False) for r in state['runs'].values()),
          '; reusing OLS:', sum(r.get('fit_reused', False) for r in state['runs'].values()))
    print('Completed 5-component task results:', 'reused unchanged' if state.get('five_component_fits_reused') else 'will reuse; no new 5-component fit')
    print('Group estimates:', state.get('component_bar_estimates', 0), '/ 450 component; ', state.get('mean_bar_estimates', 0), '/ 15 insula mean')
    print('Figures:', len(list((record/'figures').glob('*.png'))), '/ 7; magnitude review flags:', state.get('coefficient_review_flag_count', 'pending'))
    errors = [key+': '+item.get('error', 'failed') for key, item in state['runs'].items() if item['state'] == 'failed']
    if state.get('error'): errors.append(state['error'])
    errors.extend('Git publication: '+p.read_text() for p in record.glob('*publication-error.txt'))
    print('Recorded errors:', len(errors), '(up to 5 shown)')
    for error in errors[:5]: print(' ', s.short(error))
    pointer = s.ROOT/'code/insula_ica/output/task-dimensions-launcher.logpath'
    log = Path(pointer.read_text().strip()) if pointer.is_file() else None
    print('Launcher errors:', ' | '.join(s.log_errors(log)) or 'none found' if log else 'log unavailable')
    print('Results:', record)


if __name__ == '__main__': main()
