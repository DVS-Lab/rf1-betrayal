# Task responses at 3/5/7 components with insula-mean references

The requested sensitivity analysis reuses the completed unsmoothed fixed ICA
solutions. It does not rerun ICA, FEAT, preprocessing or connectivity analysis.
The completed pooled five-component task estimates are read directly from
`linux2-pooled-network-tasks-20261030T020631Z-75dba4` (or a later completed source
explicitly selected at launch). All 150 five-component group estimates must
match their saved effect/SE/CI values within 1e-12 after aggregation; otherwise
the sensitivity record fails. No five-component timecourse or task model is
refitted.

## Extraction and temporal model

- Retain the original 132 participants and exactly 490 accepted runs
  (241 Trust, 249 UGR). Fixed left/right/bilateral masks are unchanged.
- Use all existing 3/7 maps from both discovery representations (`joint_bold`
  and `native_film`). Within each mask, spatially regress the full set of
  component maps plus a spatial intercept into the exact 4D input named in the
  active FEAT `feat_files(1)` setting. Only masked voxels enter this fit; map
  voxelwise scales and existing component orientations are retained.
- The actual BOLD mean over each fixed insula mask is extracted once per run,
  independently of the component timecourses. It is a separate response, not an
  additional confound or the average of component timecourses. The same mean
  estimates are displayed alongside every decomposition for that ROI.
- All 60 new component timecourses and three ROI means are extracted during
  one new 4D read per accepted run. These new responses were not saved by the
  previous five-component task run.
- Z-score every response and full task+nuisance design column within run,
  using sample SD. Fit simultaneous OLS with the complete saved nuisance matrix,
  including the original TEDANA regressors. There is no preliminary nuisance
  regression, temporal filtering or new whitening. Constant design columns are
  recorded and zeroed by the existing fitter.
- For UGR, sum high/low pmod columns before z-scoring to fit common nonsocial
  and social slopes. Keep all four condition main effects, RT, missed and
  nuisance columns. This is the same pooled model used for the completed
  five-component results, including sub-11031 run 1.
- Trust uses the existing full design and three standardized beta contrasts:
  computer, friend and stranger reciprocate minus defect.
- Check exact saved model/confound identities, FEAT input grid/nvol/TR and
  provenance first. BOLD/FSF/TR and selected model/nuisance/mask identities must
  also match the original completed five-component task record. The corrected
  sub-10817 UGR run 1 remains TR=1.710 s.

Runs are averaged equally within participant/task, followed by equal-participant
means and pointwise 95% t CIs. Every IC has five bars: two UGR slopes and three
Trust contrasts. All 132 participants are required for every group estimate.
There are 450 component bar estimates across dimensions/bases/ROIs, plus 15
unique insula-mean bar estimates. These are descriptive summaries, without
permutation tests, significance selection or a new multiple-comparison claim.

Component numbers do not transfer between decompositions. The record includes
all spatial pairwise correlations, maximum-absolute-correlation one-to-one
matches and sign-aligned descriptive task comparisons to the five-component
basis. Matching depends only on the spatial maps. Every component, including
unmatched components, remains in the plots; no task effect selects a solution.

## Launch, resume and publication

```bash
bash code/insula_ica/start_task_dimensions.sh
```

The launcher uses the existing Linux2 `tedana-26.0.3` Python environment and
`nohup`, with eight concurrent RAM-gated 4D readers and single-threaded BLAS.
A pilot includes a Trust run, sub-11031 UGR run 1 and corrected-TR sub-10817 UGR
run 1. Expansion requires all three pilots to pass. Publication failures are
recorded without discarding completed computation or preventing cohort expansion.

Raw response arrays, standardized model/response/beta arrays, and validated
per-run extraction/fit checkpoints stay under a new ignored directory:
`derivatives/insula_ica/task-dimensions-<timestamp>/`. Original outputs remain
untouched. An interrupted OLS fit reuses its completed extraction; partial fit
files are retained with an interruption suffix. Completed fit markers validate
inputs and output hashes before reuse. Changed sources block silent reuse.
Rerunning the launcher automatically resumes the latest matching source/policy;
an exclusive lock prevents simultaneous submitters. `--resume-record` selects a
specific record and `--new-record` creates a separate analysis directory.

The Git record `code/insula_ica/records/linux2-task-dimensions-<timestamp>/`
contains accepted-run identities, map/mask copies, methods, provenance,
per-run/participant/group estimates, correspondence tables, verification of
five-component reuse, numerical review flags and seven PNGs: one plot for each
basis/dimension combination with all three ROIs, plus a mean-only plot.
Each component and mean reference receives five bars. It publishes the pilot
and final record to `insula-ica`, staging only its own record. Errors are saved
to the record and can be reviewed without dumping logs:

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/insula_ica/status_task_dimensions.py --brief
```

## Local validation (before Linux2 execution)

`test_task_dimension_sensitivity.py` checks known-source recovery at both
dimensions, spatial-intercept behavior, actual ROI mean extraction, all 63
responses fitted against an independent OLS calculation using a real synthetic
NIfTI input, both completed-fit and interrupted-fit resume without a new 4D
read, altered-cache rejection, original five-component input agreement, frozen
source protection, all 465 group summaries and preservation of five-component
summary statistics. All seven tests pass. Existing network-timecourse and
pooled-model tests are also run to check the reused estimators.

A synthetic plot smoke test generated all seven figures and 465 bar estimates;
the seven-component grid and mean-only figure were visually inspected for
labels, all five bars, CIs and legends. Synthetic plots are not scientific results
and are not included in Git. No Linux2 sensitivity result is claimed here.
