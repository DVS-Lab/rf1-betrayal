#!/usr/bin/env python3
"""Plot every primary IC and all five task bars from a completed pooled record."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import run_network_timecourse_analysis as n
from run_bold_controls import save


def plot(record):
    rows = n.g.core.read_table(record/'group_network_task_effects.tsv')
    n.g.core.require(len(rows) == 150 and all(int(r['n']) == 132 for r in rows), 'Expected all 150 N132 bar estimates')
    output = record/'overview'; output.mkdir(exist_ok=True)
    titles = ('UGR nonsocial', 'UGR social', 'Trust computer: rec - def',
              'Trust friend: rec - def', 'Trust stranger: rec - def')
    colors = ('#68a9cf', '#186da2', '#e6a658', '#bf7340', '#88532e')
    ymin = min(float(r['ci_low']) for r in rows); ymax = max(float(r['ci_high']) for r in rows)
    margin = (ymax-ymin)*.1
    for mode in n.g.inputs.MODES:
        fig, axes = plt.subplots(3, 1, figsize=(13.5, 10), sharey=True)
        for ax, roi in zip(axes, n.g.core.ROIS):
            for j, (label, title, color) in enumerate(zip(n.BAR_LABELS, titles, colors)):
                rr = [next(r for r in rows if r['mode'] == mode and r['roi'] == roi and
                           int(r['component']) == k and r['contrast'] == label) for k in range(1, 6)]
                y = np.array([float(r['effect']) for r in rr]); lo = np.array([float(r['ci_low']) for r in rr]); hi = np.array([float(r['ci_high']) for r in rr])
                ax.bar(np.arange(5)+(j-2)*.155, y, width=.145, color=color, label=title,
                       yerr=[y-lo, hi-y], capsize=2.8)
            ax.axhline(0, color='.35', lw=.9)
            ax.set(xticks=range(5), xticklabels=['IC'+str(k) for k in range(1, 6)],
                   ylabel='Standardized beta / contrast', title=roi.capitalize()+' insula',
                   ylim=(ymin-margin, ymax+margin))
            ax.spines[['top', 'right']].set_visible(False)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .952), ncol=3, frameon=False, fontsize=9)
        name = 'Joint-BOLD ICA basis' if mode == 'joint_bold' else 'Native-FILM ICA basis'
        fig.suptitle(f'{name} | all five ICs, five task bars each', fontsize=15, y=.99)
        fig.text(.5, .045, 'UGR: common high/low pmods pooled before z-scoring. Trust: partner-specific reciprocate minus defect.', ha='center', fontsize=9)
        fig.text(.5, .023, 'N=132; all 490 accepted runs retained; equal run means within participant; pointwise 95% CIs. No adjusted significance tests.', ha='center', fontsize=8.5)
        fig.text(.5, .006, 'Linear axes shared across panels. Component numbers are local to each ROI solution.', ha='center', fontsize=8)
        fig.tight_layout(rect=(0, .07, 1, .90)); fig.savefig(output/f'{mode}_all_ICs_five_bars.png', dpi=170); plt.close(fig)
    save(output/'provenance.json', dict(source=n.g.core.identity(record/'group_network_task_effects.tsv'),
        plotted_values=150, bars_per_component=5, components_per_ROI=5, ROIs=list(n.g.core.ROIS),
        bases=list(n.g.inputs.MODES), scale='Linear and shared; no clipping, outlier removal or refitting',
        titles=titles))
    print('Complete five-bar overviews:', output)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--record', type=Path); args = p.parse_args()
    plot(args.record or sorted(n.RECORDS.glob('linux2-pooled-network-tasks-*'))[-1])
