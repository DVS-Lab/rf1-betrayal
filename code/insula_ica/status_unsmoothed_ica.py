#!/usr/bin/env python3
"""Bounded corrected unsmoothed ICA progress and errors; no imaging reads."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import status_insula_streams as s


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--brief',action='store_true');p.parse_args()
    records=s.ROOT/'code/insula_ica/records'
    choices=sorted(p for p in records.glob('linux2-unsmoothed-group-ica-*') if (p/'state.json').is_file())
    if not choices:print('No corrected unsmoothed ICA record found.');return
    record=choices[-1];state=s.read(record/'state.json');rows=s.processes()
    runners,active,counts=s.pipeline(rows or [],('run_unsmoothed_ica.py',),Path(state['output_root']))
    errors=[]
    for group in ('runs','fits'):
        errors.extend(k+': '+v.get('error','failed') for k,v in state.get(group,{}).items() if v.get('state')=='failed')
    if state.get('error'):errors.append(state['error'])
    if (record/'publication-error.txt').is_file():errors.append('Git publication: '+(record/'publication-error.txt').read_text())
    fits=state.get('fits',{});done=[k for k,v in fits.items() if v.get('state')=='complete'];auto=sum(k.endswith('/dim-0') for k in done)
    print('Record:',record.name)
    print('Runner:','unknown' if rows is None else ', '.join(f'PID {r["pid"]}, elapsed {r["elapsed"]}' for r in runners) if runners else 'not running')
    print('State:',state['state'])
    print('Paired residual runs prepared:',sum(v.get('state')=='complete' for v in state.get('runs',{}).values()),'/490')
    print(f'ICA fits complete: {len(done)}/24 (fixed {len(done)-auto}/18; automatic {auto}/6)')
    if rows is not None:print(f'Active MELODIC {counts["melodic"]}; FEAT {counts["feat"]}; FILM {counts["film_gls"]}')
    for key in done:
        if key.endswith('/dim-0'):print('Automatic components:',key,'=',fits[key]['actual_components'])
    print('Recorded errors:',len(errors),'(up to 6 shown)')
    for error in errors[:6]:print(' ',s.brief_error(error))
    log=s.launcher('unsmoothed-group-ica','rf1-unsmoothed-group-ica.*');tail=s.log_errors(log)
    print('Launcher errors:',' | '.join(tail) if tail else 'none found' if log else 'log unavailable')
    saved=s.ROOT/'code/insula_ica/output'/('unsmoothed-ica-status-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));saved.mkdir(parents=True)
    (saved/'snapshot.json').write_text(json.dumps(dict(record=str(record),state=state,errors=errors,runner_processes=runners,
        active_processes=active,launcher_log=str(log) if log else None,launcher_errors=tail),indent=2)+'\n')
    print('Full snapshot:',saved/'snapshot.json')


if __name__=='__main__':main()
