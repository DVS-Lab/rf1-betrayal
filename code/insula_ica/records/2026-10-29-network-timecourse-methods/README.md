# Direct network timecourse task analysis

The previous spatial projection of voxelwise task estimates is withdrawn. Its
code, six published task artifacts, numerical findings and directional task
conclusions were deleted from the working branch on 2026-10-29. Git history is
retained; no history rewrite or change to original FEAT statistics is intended.
Completed ICA fits, anatomical descriptions and descriptive background dual
regression are retained. No previous IC3/IC4 task direction is carried forward.

## Paper and framework

Dobryakova & Smith (2022), *Reward enhances connectivity between the ventral
striatum and the default mode network*, NeuroImage 258:119398,
[doi:10.1016/j.neuroimage.2022.119398](https://doi.org/10.1016/j.neuroimage.2022.119398).
Sections 2.5 and 3.1, and Figure 1, distinguish spatially extracted network
timecourses, temporal modeling of their task responses, and nPPI interaction
models for task-dependent connectivity. This paper uses whole-brain reference
maps; our spatial model uses the completed restricted-insula ICA maps.

The lab's [Leech-style script](https://github.com/DVS-Lab/duke-archive/blob/master/leech12_dual_regression)
extracts all component timecourses inside the restricted mask and temporally
regresses them into whole-brain data. Its
[confound extension](https://github.com/DVS-Lab/duke-archive/blob/master/leech12_dual_regression_conf)
appends nuisance regressors at the temporal stage. Standard
[FSL dual regression](https://fsl.fmrib.ox.ac.uk/fsl/docs/resting_state/dualregression.html)
uses spatial maps -> individual timecourses -> individual connectivity maps.
The task-response model uses the extracted timecourses as responses to the task
design; connectivity uses them as predictors of voxel timecourses. These are
distinct models and outputs.

## Replacement task estimator

1. Freeze all six completed primary five-component solutions: joint-BOLD/native-
   FILM discovery bases by left/right/bilateral insula. Retain every component;
   component indices are local to each solution. Keep the existing sign rule
   (largest absolute spatial loading positive). No effect determines selection
   or orientation. Discovery remains based on the already agreed residuals;
   no new ICA is estimated here.
2. Read each run's corrected, unsmoothed, **task-preserved** `filtered_func_data`.
   Spatially regress all five maps together, with a spatial intercept, inside
   their defining mask to estimate T x 5 network timecourses. Centering the
   timecourses removes only a constant baseline. Do not extract task responses
   from `res4d` or the original-scale full-model residuals.
3. Fit each network response to the **complete saved task+nuisance design** in a
   single temporal GLS model. This includes the exact preserved TEDANA/confound
   matrix. Source SHA, numeric nuisance equality, completed source QC, grid,
   timing, model/contrast coding and estimability are checked. The corrected
   UGR 1.710-second model is retained; other runs remain 1.615 seconds. No new
   high-pass, nuisance removal, HRF convolution or sequential regression occurs.
4. Use FILM's own component-level prewhitening, without spatial autocorrelation
   smoothing. Five responses occupy five independent artificial voxels in one
   small input, so each receives its own autocorrelation estimate. `--sa` and
   `--noest` are absent. A permissive `--thr` prevents dropping mean-zero
   timecourses. The new contrast outputs come from this temporal model; no
   existing voxelwise task map is read by the task estimator. Save unscaled
   stage-1 timecourses, fitted contrast variances, commands and a direct temporal
   OLS diagnostic. Do not normalize task responses by residual SD or claim PSC.
5. For Trust, retain the original summed C/F/S contrast units and report
   betrayal minus reciprocity. Estimate the betrayal and reciprocity condition
   means separately using their actual model weights (C/F/S average). Their
   difference equals one third of the primary summed contrast, not the entire
   contrast. UGR uses its stored offer-pmod contrast without an inferred fairness
   sign. All nuisance coefficients have zero task-contrast weights.
6. Average runs equally within participant/task, then test participants, N=132.
   Retain all 60 primary tests across both map bases, three ROIs, five components
   and two tasks, with Bonferroni-60. Condition bars show pointwise between-person
   95% CIs and have no separate significance claims. All primary estimates and
   nulls are published. This is exploratory same-cohort map discovery, not an
   independent replication or a formal IC3-vs-IC4/task/hemisphere interaction.

**Both discovery bases are applied to the same unwhitened task-preserved BOLD.**
`native_film` labels where an ICA spatial basis was discovered, not the data
representation entering this new temporal response model. Voxelwise original
FILM whitening is not reused for these extracted timecourses.

## Connectivity and scope

Existing connectivity estimates are background dual-regression results: all
five spatially derived residual timecourses enter a simultaneous whole-brain
regression on the same residual representation. The original joint FEAT fit
already removed task and nuisance effects. They remain separate from this new
task-preserved stage-1/timecourse analysis. No new PPI or claim of condition-
dependent connectivity is introduced. The attached paper's nPPI extension
would require network x task interactions and an explicit whole-brain temporal
model; a main task effect on a network is not such an interaction.

## Execution and outputs

```bash
bash code/insula_ica/start_network_timecourses.sh
```

Detached from the terminal. Runs Trust, ordinary UGR and the corrected-TR UGR
pilot first, then automatically expands with eight RAM-gated 4D readers. This
does no new ICA, fMRIPrep, FEAT prestats or voxelwise model fitting. It performs
2940 small network-level GLS fits (490 runs x 6 map bases), using existing
complete inputs. Partial fits are preserved and never silently overwritten.

```bash
python code/insula_ica/status_network_timecourses.py --brief
# For a failed/interrupted record, retain completed runs and retry failed ones:
python code/insula_ica/run_network_timecourse_analysis.py \
  --resume-record /absolute/path/to/linux2-network-timecourses-RECORD \
  --retry-failed --jobs 8 --publish
```

Linux2 Python is `/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python`.
Small group/participant/run result tables, plots and provenance publish under
`records/linux2-network-timecourses-*`. Extracted timecourses and model outputs
stay under a separate ignored `derivatives/insula_ica/network-timecourses-*`
root. The launcher removes the specifically withdrawn ignored exploration task
outputs and strips their estimates from completion markers, while retaining all
FEAT maps, ICA fits, connectivity maps and 4D data. Git pull removes the tracked
withdrawn artifacts. The cleanup log is retained under ignored `output/`.

## Validation completed locally

Six replacement tests passed: safe withdrawal cleanup preserving source FEAT and
connectivity maps; task-signal recovery with correlated nuisance; actual
condition-contrast construction and EV-order validation; participant/run
aggregation and all 60 primary tests; real-FSL GLS estimation of five independent
component timecourses against a single-timecourse reference; and a real-FSL
three-run 4D pilot/resume exercising all six map bases without source task maps.
A separate real-FSL stage-1 check matched `fsl_glm --demean` (pooled r > 0.9999999999; maximum absolute float32 difference 0.000033).
Five background-connectivity tests and nine retained inventory/numerical tests
also passed. Shell syntax and whitespace checks passed. The bar-plot layout was
rendered with synthetic values in temporary storage only. No synthetic results
are published as participant findings. Linux2 RF1 pilot/final results remain
pending.
