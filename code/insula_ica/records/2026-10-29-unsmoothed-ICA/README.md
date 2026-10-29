# Corrected unsmoothed paired-scale insula ICA

The timing-corrected sm-0 registry completed 490/490 FEATs (132 participants,
Trust 241 and UGR 249 runs). This launcher uses that registry and requires the
authorized sub-10817 UGR case's regenerated 1.710-second task design. All source
completion outputs and model/nuisance QC are verified before preparation. The
historical model path is retained as provenance, while each validated sm-0 full
model supplies the parameter estimates and native residuals for this analysis.

```bash
bash code/insula_ica/start_unsmoothed_ica.sh
```

There are 24 fits: left/right/bilateral × joint-BOLD/native-FILM residuals ×
3/5/7/automatic dimensions. Five remains primary. Native MELODIC variance
normalization is enabled consistently in this stream, alongside the documented
runwise SD/energy balancing. Each saved input has a uniform ×1000 multiplier to
avoid the absolute SD floor in classic MELODIC normalization; its mean is not
added back. This normalization policy differs from the older manual-only baseline
and must be labeled when comparing those exploratory results.

Each ROI gets its own physically masked 4D files, verified exactly zero outside
that ROI. The common support is recalculated across all sources and both scales.
All accepted frames remain in participant/task blocks; run contributions are
balanced as before. The corrected participant has one accepted UGR run at 1.710
seconds, so its block retains that header. Other blocks retain 1.615 seconds.
Different TRs within a block cause an explicit failure instead of interpolation
or relabeling. Input signatures and contribution manifests retain the TRs.

This is spatial temporal-concatenation ICA. No extra temporal filtering/regression
or resampling is introduced. A global `--tr=1.615` override and MELODIC web report
are omitted. Review figures contain spatial maps only; mixed-TR concatenated
time/frequency plots are not interpreted. FSL documents TR's role in plot units
on its [MELODIC page](https://fsl.fmrib.ox.ac.uk/fsl/docs/resting_state/melodic.html).
Later time-domain inference must continue using each run's actual acquisition TR.

The launcher queues 64 residual-preparation workers and 18 fixed ICA fits, bounded
by RAM reservations. Six automatic fits run serially, using classic PCA, eight
threads and a 32 GiB reservation plus eight GiB headroom. This prevents six
simultaneous large classic-PCA allocations; it is not a claim that the prior
segmentation fault's cause was established. Auto fits run even if a fixed fit
fails. The same five deterministic starts, 10,000 iterations, pow3 objective and
epsilon 0.00005 are used. Neither objective nor convergence tolerance changes on
failure. Successful outputs and bounded retry logs are preserved.

Native residuals and PEs already exist: this launcher requests no fMRIPrep, FEAT,
FILM, TEDANA, or Warpkit. It writes arrays, masked inputs and fits into a dedicated
`derivatives/insula_ica/unsmoothed-group-ica-*` root, preserving all old analysis
outputs. Running it again reuses the latest unsmoothed root and completed verified
preparation/fits. Unknown partial ICA outputs require inspection rather than
automatic deletion/replacement. A lock prevents duplicate submitters.

Records publish to `insula-ica` at completion/failure. Full success generates all
spatial figures, residual-scale comparisons and dimensionality stability. Partial
success still produces spatial review figures for completed fits and records all
failures. Whole-brain connectivity, task inference, task-preserved OLS sensitivity,
and the manuscript's corrected group/PPI propagation remain separate steps.

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/insula_ica/status_unsmoothed_ica.py --brief
```

## Validation

31 selected tests passed: 25 existing group/resume/task-preserved regressions and
six new tests. They cover unchanged legacy commands, all 24 fit orchestration,
continued automatic fits after a fixed-fit failure, the corrected 490-run source
registry and rejection of missing corrected-design QC, physical ROI masking for
all 1584 input files, preserved 1.615/1.710 headers, and refusal to combine differing
TRs within one block. The shell launcher passes syntax checks.

A real local MELODIC probe with 264 physically masked synthetic inputs, native VN,
both TRs, no common-TR override and no web report converged to three components.
Its output is `/private/tmp/rf1-mixed-TR-ICA-6cvx133x/probe-result.json` (local FSL
6.0.7.19; linux2 records FSL 6.0.7.18). This is software validation, not completion
of the RF1 fits or verification of full-data automatic-PCA memory requirements.
