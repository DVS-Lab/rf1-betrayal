# Priority task review: standardized OLS and five bars per subregion

Current policy: `FEAT-input-network-timecourse-zscored-OLS-v2`, clarified by the
user on 2026-10-29. This replaces the earlier unstandardized GLS proposal. The
original spatial task-map projections remain deleted and withdrawn.

## Task branch (run first)

1. Freeze all six completed primary five-component ICA bases: left/right/bilateral
   insula by joint-BOLD/native-FILM discovery representation. All components and
   existing sign orientations are retained. No new ICA is estimated.
2. Read the exact 4D `feat_files(1)` input in each run's active saved `design.fsf`.
   Verify its grid, volumes and run-specific TR against the actual L1 design.
   This is the input to FEAT, **not** `filtered_func_data`, `res4d` or the full-
   model original-scale residuals. Restrict the data matrix to voxels in the
   defining ROI before spatial regression. No outside-mask voxel enters the fit.
3. Simultaneously spatially regress all five component maps with a spatial
   intercept into that masked input. The resulting network timecourses retain
   task and nuisance effects. Only a temporal constant baseline is removed.
4. Columnwise z-score each network timecourse and **every column** of the saved
   complete L1 task+nuisance design, within run, using sample SD (`ddof=1`).
   Constant design columns are represented by zero and explicitly listed. No
   meaningful task or nuisance regressor is excluded. Verify source confound SHA
   and numeric equality with the exact historical saved TEDANA/confound matrix
   before standardization. Preserve the corrected UGR 1.710-second design and
   the other runs' 1.615-second timing.
5. Fit `Z(TC) = Z(X_full) B + error` with OLS (`pinv(Z(X_full)) @ Z(TC)`). No
   temporal whitening, FILM, additional filtering or sequential task/nuisance
   cleaning. An intercept is zero because all inputs are centered. Check that
   requested task contrasts are estimable. Save raw extracted timecourses, both
   standardized matrices, the full standardized coefficient matrix, means/SDs,
   source identities and numerical rank.

## Exactly five plotted values per subregion

| Bar | Operation on fitted standardized coefficients |
| --- | --- |
| UGR nonsocial | `(beta_nonsocial_high_pmod + beta_nonsocial_low_pmod) / 2` |
| UGR social | `(beta_social_high_pmod + beta_social_low_pmod) / 2` |
| Trust C | `beta_C_rec - beta_C_def` |
| Trust F | `beta_F_rec - beta_F_def` |
| Trust S | `beta_S_rec - beta_S_def` |

Standardize individual EV columns **before** fitting and forming these averages
or contrasts. Do not average the high/low design columns before standardizing.
Keep all four individual UGR pmod estimates in the result tables. Pmod direction
follows the original L1 coding; no inferred unfairness sign flip is applied.

Calculate the five values within run, average runs equally within participant
and task, then plot equal-weight group means over N=132 with pointwise 95%
between-participant t CIs. In particular, CIs for an average or difference use
that participant's actual combined estimate, retaining the paired covariance.
There are 150 plotted group values: two discovery bases x three ROIs x five
components x five bars. The component indices are local to each basis/ROI.

Produce one grouped five-bar figure for each ROI/basis (six figures) plus
bilateral IC3/IC4 close-ups (two figures). This is the priority descriptive
review; no significance stars, permutation tests or new multiple-testing family
are introduced here. Standardization does not by itself establish a formal
between-task or between-component difference. Discovery uses this same cohort,
so results remain exploratory.

## Background connectivity branch (separate; group inference later)

Spatial regression on **masked residualized** data yields background component
timecourses. Simultaneously regress all component timecourses into the matching
whole-brain residualized data to obtain run-level connectivity maps, separately
for the original-BOLD and native FILM residual representations. The existing
full FEAT fit already removed task and nuisance effects. These timecourses are
not used for task-response estimation. Average runs within participant/task
before treating participants as independent observations.

`randomise` is for group inference on the participant connectivity maps; it does
not perform the run-level timecourse-to-whole-brain regression. The task-bar
launcher does not launch either connectivity regression or `randomise`.

## Relation to the supplied paper

Dobryakova & Smith (2022), *Reward enhances connectivity between the ventral
striatum and the default mode network*, NeuroImage 258:119398,
[doi:10.1016/j.neuroimage.2022.119398](https://doi.org/10.1016/j.neuroimage.2022.119398).
Sections 2.5/3.1 and Figure 1 distinguish spatially extracted network timecourses,
temporal task-response modeling, and nPPI models for task-dependent connectivity.
The user's current instruction specifies standardized OLS for our task branch.
No nPPI interactions are introduced at this stage. The lab's
[restricted dual-regression script](https://github.com/DVS-Lab/duke-archive/blob/master/leech12_dual_regression)
and [confound extension](https://github.com/DVS-Lab/duke-archive/blob/master/leech12_dual_regression_conf)
provide the corresponding connectivity framework.

## Run on Linux2

```bash
bash code/insula_ica/start_network_timecourses.sh
```

Detached automatically. Tests one Trust, one standard-TR UGR and the corrected-
TR UGR run, then expands to all 490 accepted runs using eight RAM-gated readers.
The computation is now just spatial regression plus small OLS fits, with no
external statistical subprocess per component or run. Pilot/final/failure
records publish to GitHub. This makes new task outputs and retains all existing
ICA and background-connectivity work. An older GLS-policy record cannot be
resumed under this new estimator; start a fresh task record instead.

```bash
python code/insula_ica/status_network_timecourses.py --brief
# Resume only a record with the current standardized-OLS policy:
python code/insula_ica/run_network_timecourse_analysis.py \
  --resume-record /absolute/path/to/linux2-network-timecourses-RECORD \
  --retry-failed --jobs 8 --publish
```

Linux2 Python: `/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python`.
Small tables, eight task plots and provenance publish under
`records/linux2-network-timecourses-*`. Raw timecourses, standardized matrices
and full fitted coefficients remain in the separate ignored
`derivatives/insula_ica/network-timecourses-*` root. Failed partial outputs are
preserved rather than overwritten. Completed OLS runs can be resumed without
re-extraction if all frozen input/output identities still match.

## Validation

Six targeted tests passed: safe withdrawal cleanup; spatial extraction with
correlated task/nuisance signal; all five bar weights and EV ordering; actual
columnwise standardization and joint OLS versus an independent least-squares
solution, including constant columns; run-before-participant aggregation and
all 150 planned bar values; and a three-run 4D pilot/resume that proves the exact
FEAT input is used despite a deliberately different `filtered_func_data`.
Linux2 participant results and real-data bars are pending execution.
