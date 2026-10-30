#!/usr/bin/env python3
"""Audit cached standardized task designs, without 4D reads, exclusion or refits."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import uuid

import numpy as np

import run_network_timecourse_analysis as n
from run_bold_controls import save, publish


def design_diagnostics(zx, task_columns):
    active = np.linalg.norm(zx, axis=0) > 1e-10
    values = zx[:, active]; singular = np.linalg.svd(values, compute_uv=False)
    diagnostic = dict(rows=len(zx), columns=zx.shape[1], nonconstant_columns=int(active.sum()),
        singular_values=singular.tolist(), condition_nonconstant=float(singular[0]/max(singular[-1], 1e-300)),
        rank_at_relative_cutoffs={str(c):int((singular > singular[0]*c).sum()) for c in (1e-10, 1e-8, 1e-6)})
    columns = []
    for j in range(task_columns):
        target = zx[:, j]; others = np.delete(zx, j, axis=1)
        norm = np.linalg.norm(target)
        residual = target-others@np.linalg.pinv(others, rcond=1e-10)@target
        columns.append(dict(index=j, conditional_energy_fraction=float((residual@residual)/(norm*norm)) if norm > 0 else None,
                            constant=bool(norm == 0)))
    diagnostic['task_columns'] = columns
    return diagnostic


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--record', type=Path)
    parser.add_argument('--publish', action='store_true'); args = parser.parse_args()
    source = args.record or sorted(n.RECORDS.glob('linux2-network-timecourses-*'))[-1]
    state = json.loads((source/'state.json').read_text()); output = Path(state['output_root'])
    rows = json.loads((source/'accepted-inputs.json').read_text())
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    record = n.RECORDS/('linux2-network-design-audit-'+stamp); record.mkdir()
    summaries = []; details = []
    for row in rows:
        feat = Path(row['active_feat']); settings = n.g.core.fsl_settings(feat/'design.fsf')
        folder = output/row['key']/'joint_bold_bilateral'
        zx = np.load(folder/'standardized_design.npy'); diag = design_diagnostics(zx, int(settings['fmri(evs_real)']))
        weakest = min((r['conditional_energy_fraction'] for r in diag['task_columns'] if r['conditional_energy_fraction'] is not None), default=1.)
        summaries.append(dict(key=row['key'], condition=diag['condition_nonconstant'], min_task_conditional_energy=weakest))
        if weakest < 1e-8 or diag['condition_nonconstant'] > 1e6 or row['subject'] == '11031':
            detail = dict(key=row['key'], diagnostic=diag, design=n.g.core.identity(feat/'design.mat'),
                          standardization=json.loads((folder/'standardization.json').read_text()),
                          task_EV_titles=[settings.get(f'fmri(evtitle{k})', '') for k in range(1, int(settings['fmri(evs_real)'])+1)])
            events = []
            for k in range(1, int(settings['fmri(evs_orig)'])+1):
                path = Path(settings.get(f'fmri(custom{k})', ''))
                if not path.is_file(): continue
                values = np.loadtxt(path, ndmin=2)
                if values.shape[1] < 3: continue
                events.append(dict(EV=k, title=settings.get(f'fmri(evtitle{k})'), path=str(path), sha256=n.g.core.sha(path),
                    events=len(values), amplitude_min=float(values[:, 2].min()), amplitude_max=float(values[:, 2].max()),
                    amplitude_SD=float(values[:, 2].std()), unique_amplitudes=int(len(np.unique(values[:, 2])))))
            detail['event_amplitude_summaries'] = events
            if row['subject'] == '11031':
                sensitivities = []
                for mode in n.g.inputs.MODES:
                    zy = np.load(output/row['key']/f'{mode}_bilateral/standardized_timecourses.npy')
                    cs, labels = n.contrasts_for(row, settings, zx, n.g.core.fsl_matrix(feat/'design.con'))
                    for cutoff in (1e-10, 1e-8, 1e-6):
                        beta = np.linalg.pinv(zx, rcond=cutoff)@zy
                        sensitivities.append(dict(mode=mode, cutoff=cutoff, labels=labels, diagnostic_estimates=(cs@beta).tolist(),
                                                  note='Diagnostics only; no result replacement or exclusion.'))
                detail['cutoff_diagnostics'] = sensitivities
            details.append(detail)
    n.g.core.write_table(record/'design_summary.tsv', summaries); save(record/'flagged_designs.json', details)
    save(record/'state.json', dict(state='design_diagnostic_review_ready', source_record=str(source), runs=len(rows), detailed_runs=len(details),
        no_4D_reads=True, no_exclusions=True, no_refitting_of_production_results=True))
    print('Designs audited:', len(rows), '| detailed flags/control runs:', len(details), '| record:', record)
    if args.publish: publish(record, 'Publish standardized task-design diagnostics without exclusions or production refits')


if __name__ == '__main__': main()
