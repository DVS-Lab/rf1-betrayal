#!/usr/bin/env python3
"""Compact matched 5-mm arm status, with no participant image reads."""
import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code/insula_ica'))
import status_insula_streams as s


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--brief',action='store_true');p.parse_args()
    choices=sorted(p for p in (ROOT/'code/functional_cluster_ica/records').glob('linux2-functional-sm5-*') if (p/'state.json').is_file())
    if not choices:print('No matched sm-5 record yet.');return
    record=choices[-1];state=s.read(record/'state.json');rows=s.processes()
    runners,_,counts=s.pipeline(rows or [],('code/functional_cluster_ica/smoothed.py',),Path(state['output_root']))
    print('Record:',record.name)
    print('Runner:','unknown' if rows is None else ', '.join(f'PID {r["pid"]}, elapsed {r["elapsed"]}' for r in runners) if runners else 'not running')
    print('State:',state['state'])
    for key,title in [('feat_runs','Matched sm-5 FEAT + checks'),('residual_runs','Paired residual runs'),('component_runs','Functional IC task fits')]:
        values=state.get(key,{})
        print(title+':',sum(v['state']=='complete' for v in values.values()),'/490; failed:',sum(v['state']=='failed' for v in values.values()))
    print('Verified sm-5 FEATs reused:',sum(v.get('reused_completed_feat',False) for v in state.get('feat_runs',{}).values()))
    done=[v for v in state.get('fits',{}).values() if v['state']=='complete']
    print('ICA complete:',len(done),'/8; fixed:',sum(v['dimensionality']!=0 for v in done),'/6; automatic:',sum(v['dimensionality']==0 for v in done),'/2')
    print('Unsmoothed ICA reruns:',state.get('existing_unsmoothed_ICA_reruns',0),'| mean refits:',state.get('mean_fits_rerun',0))
    print('Active processes:',', '.join(k+' '+str(v) for k,v in counts.items()))
    errors=[]
    for section in ('feat_runs','residual_runs','fits','component_runs'):
        errors.extend(k+': '+v.get('error','failed') for k,v in state.get(section,{}).items() if v['state']=='failed')
    if state.get('error'):errors.append(state['error'])
    if (record/'publication-error.txt').is_file():errors.append('Git publication: '+(record/'publication-error.txt').read_text())
    print('Recorded errors:',len(errors),'(up to 6 shown)')
    for error in errors[:6]:print(' ',s.short(error))
    pointer=ROOT/'code/insula_ica/output/functional-sm5-launcher.logpath';log=Path(pointer.read_text().strip()) if pointer.is_file() else None
    print('Launcher errors:',' | '.join(s.log_errors(log)) or 'none found' if log else 'log unavailable')
    print('Results:',record)


if __name__=='__main__':main()
