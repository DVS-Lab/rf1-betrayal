#!/usr/bin/env python3
"""Review published task tables, retain every participant, and show Trust clearly."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import run_network_timecourse_analysis as n
from run_bold_controls import save


def coefficient_flags(rows, threshold=5.):
    # A review trigger, not an exclusion rule or a bound on standardized beta.
    # Large partial coefficients can occur with highly collinear predictors.
    return [r for r in rows if abs(float(r['effect'])) > threshold]


def review(record):
    rows = n.g.core.read_table(record/'group_network_task_effects.tsv')
    runs = n.g.core.read_table(record/'run_network_task_effects.tsv')
    flags = coefficient_flags(runs); output = record/'review'; output.mkdir(exist_ok=True)
    save(output/'coefficient-review.json', dict(threshold=5., flagged_run_effects=flags,
        run_keys=sorted({f'sub-{r["subject"]}_task-{r["task"]}_run-{r["run"]}' for r in flags}),
        status='UGR social bars held for design/numerical review' if flags else 'No magnitude review flags',
        exclusions_applied=False, source_group_table=n.g.core.identity(record/'group_network_task_effects.tsv'),
        meaning='Magnitude flags request investigation; standardized partial coefficients are not bounded by one. No clipping, deletion, refitting, or outlier exclusion.'))
    colors = ('#ddab62', '#ba713f', '#754c32'); labels = ('C: rec - def', 'F: rec - def', 'S: rec - def')
    def draw(ax, mode, roi, components):
        for j, (contrast, label, color) in enumerate(zip(n.BAR_LABELS[2:], labels, colors)):
            rr = [next(r for r in rows if (r['mode'], r['roi'], int(r['component']), r['contrast']) == (mode, roi, k, contrast)) for k in components]
            y = np.array([float(r['effect']) for r in rr]); lo = np.array([float(r['ci_low']) for r in rr]); hi = np.array([float(r['ci_high']) for r in rr])
            ax.bar(np.arange(len(components))+(j-1)*.23, y, width=.21, color=color, label=label,
                   yerr=[y-lo, hi-y], capsize=4)
        ax.axhline(0, color='.35', lw=1)
        ax.set(xticks=range(len(components)), xticklabels=['IC'+str(k) for k in components],
               ylabel='Standardized beta: reciprocate - defect')
        ax.spines[['top', 'right']].set_visible(False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, mode, title in zip(axes, n.g.inputs.MODES, ('Joint-BOLD ICA basis', 'Native-FILM ICA basis')):
        draw(ax, mode, 'bilateral', (3, 4)); ax.set_title(title)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle('Trust responses | bilateral insula IC3 and IC4', fontsize=15)
    fig.text(.5, .045, 'N=132; all participants retained; 95% pointwise participant CIs; no multiplicity-adjusted significance claim.', ha='center', fontsize=9)
    fig.text(.5, .012, 'Positive: reciprocity > defection. Negative: defection > reciprocity. UGR social estimates require separate numerical review.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .1, 1, .94)); fig.savefig(output/'bilateral_IC3_IC4_Trust.png', dpi=180); plt.close(fig)
    fig, axes = plt.subplots(3, 2, figsize=(12, 11), sharey=True)
    for i, roi in enumerate(n.g.core.ROIS):
        for j, mode in enumerate(n.g.inputs.MODES):
            draw(axes[i, j], mode, roi, range(1, 6)); axes[i, j].set_title(f'{roi} | {mode}')
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle('Trust | all primary insula components | direct standardized OLS', fontsize=14)
    fig.text(.5, .015, 'N=132; equal run means within participant; pointwise 95% CIs. All estimates retained; no outlier exclusions or significance stars.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .96)); fig.savefig(output/'all_ROIs_Trust.png', dpi=150); plt.close(fig)
    return flags


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--record', type=Path)
    args = p.parse_args()
    record = args.record or sorted(n.RECORDS.glob('linux2-network-timecourses-*'))[-1]
    flags = review(record)
    print('Review:', record/'review', '| magnitude review flags:', len(flags))
