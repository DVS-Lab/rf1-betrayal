# Functional-cluster mean and exploratory masked ICA

This stream uses the documented current **127-voxel**
`masks/seed-AIns-clusterthresh.nii.gz`, SHA256
`66be24139534438a8f2efd49d071afcfa21010df5a578b955917a22b71d234d0`.
It is the surviving functional cluster candidate, not a proven binary copy of
the historical PPI seed. It lies mainly in left frontal orbital/inferior frontal
and opercular cortex; its anatomical-insula overlap is 2/127 voxels under the
Harvard–Oxford max-probability 25% definition. The existing anatomical-insula
analyses remain a different ROI analysis.

## Linux2

```bash
cd /ZPOOL/data/projects/rf1-betrayal
git pull --ff-only origin insula-ica
bash code/functional_cluster_ica/start.sh
```

The launcher selects the existing `tedana-26.0.3` Python, initializes the usual
FSL environment if necessary, and detaches with `nohup`. It starts 32 RAM-gated
extraction/temporal workers, six fixed ICA slots (four threads each), and runs
the two automatic-dimension fits serially (eight threads and 32 GiB reservation,
plus 8 GiB headroom). Its lock prevents duplicate active analyses. Repeating
the same command resumes matching completed checkpoints automatically.

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/functional_cluster_ica/status.py --brief
```

## Execution order

1. Require the latest completed `linux2-task-dimensions-*` record, exact N132
   identities/490 accepted runs, corrected sm-0 registry, and matching saved
   task models, nuisance files and raw FEAT inputs. Pilot first on a Trust run,
   sub-11031 UGR run 1 and sub-10817 UGR run 1. Keep the established pooled
   high/low UGR pmod model and corrected run-specific TR.
2. Read the exact task-preserved `feat_files(1)` BOLD once per new run and cache
   **only the original 127 ROI voxels**. Fit the actual ROI mean to the complete
   z-scored task+nuisance design. No prior nuisance or task subtraction is used
   for these task-response fits. Every series gets five bars: pooled UGR
   nonsocial/social fairness slopes, then Trust C/F/S reciprocate-minus-defect.
3. Reuse all **465 existing anatomical-insula group estimates unchanged**.
   Produce six 3/5/7-component figures (two residual bases, all L/R/bilateral
   panels), each containing the functional mean, anatomical mean and every IC.
   Also produce a four-ROI mean-comparison plot. Publish these results first.
4. Verify retained sm-0 filtered data/design/PE/native residual identities
   against the completed anatomical ICA's original input snapshots. Extract
   two functional residual representations: centered filtered BOLD minus the
   saved full design × existing FILM PEs, and the existing native prewhitened
   `res4d`. This code has **no FILM reconstruction path**. It adds no regression,
   whitening, filter or smoothing. Missing sources stop this stage.
5. Record common fitted/finite/nonconstant support across all runs and both
   scales; retain all subjects/runs. ICA uses that fixed common subset of the
   127-voxel candidate, with the existing minimum 80% support gate. Report exact
   retained voxel count; do not label an eroded subset as the original 127.
   The mean remains an average over the original 127 raw-BOLD voxels.
6. Build 264 participant/task blocks per residual scale. Demean/SD-scale each
   run and balance run energy as in the anatomical pipeline; retain all frames
   and acquisition TRs without interpolation. **Physically zero all 4D values
   outside the specific functional mask**, crop the array without changing
   coordinates, and verify the written files before MELODIC.
7. Run eight spatial temporal-concatenation fits: two residual representations
   × **3/5/7/automatic (`-d 0`)**. Use the existing deterministic five-start,
   10,000-iteration retry helper, `-a concat`, native MELODIC variance
   normalization enabled, `--disableMigp`, and no common-TR override. Full
   argument arrays and attempts are saved. Failed fits do not discard completed
   ones or the mean plots. Automatic fits are attempted even after a fixed-fit
   failure, and do not run alongside one another.
8. For each completed functional ICA, jointly regress all component maps plus
   a spatial intercept onto **masked task-preserved FEAT-input BOLD**. Use the
   cached 127-voxel values; no extra 4D reads. Regress resulting timecourses
   onto the same complete standardized temporal model. Save all components,
   spatial-map figures, five-bar plots with functional/left-insula mean
   references, and 3/7-vs-5 spatial correspondence tables. Automatic solutions
   also get paginated task/spatial figures. No COPE projection is used.

Runs are averaged within participant/task, then participant means receive equal
weight. Error bars are descriptive, pointwise 95% participant CIs. The functional
ROI was selected from task effects, so these task associations are exploratory
characterization, not independent confirmation of the effects defining the ROI.
ICs are signed spatial loadings, not disjoint parcels; labels are local to each
solution. No whole-brain connectivity, randomise, PPI, behavioral selection or
new significance testing is included in this stage.

## Storage, publication and reuse

Small records, masks/component maps, tables, plots, commands and failures publish
to `code/functional_cluster_ica/records/linux2-functional-cluster-*` on GitHub
branch `insula-ica`. Only the current record is staged by the publisher. All
participant 4D images, masked arrays, temporal caches, MELODIC outputs and session
blocks stay in a separate ignored
`derivatives/insula_ica/functional-cluster-*` root. No files are written beside
original FEAT outputs and no original designs or masks are modified.

Mean plots publish before the ICA stage; final/partial/failure records publish
at exit. A publication failure is retained in `publication-error.txt`; the status
script reports it. Repeating the launcher retries publication and resumes safely.
Source/mask/code provenance and output identities must match. Completed raw ROI
caches, mean fits, residual arrays, physical blocks, ICA solutions and each
basis's temporal fits are reused. Partial files are retained with an interrupted
suffix rather than removed. Changed inputs are rejected rather than mixed.

## Validation

Six local tests cover source/content mismatch rejection; exact raw/fit cache
reuse; simultaneous spatial recovery with a spatial intercept at 3/5/7/11
components; standardized full-design OLS; retained-PE residual construction with
FILM recovery forbidden; native-TR physical masking of all 528 synthetic
participant/task × residual-mode blocks and their reuse; and all seven new mean
figures using the actual 465 saved anatomical estimates without mutating them.
These tests do not assert that a real Linux2 cohort fit or new ICA has already
run. The production pilots and MELODIC convergence checks remain part of the
Linux2 execution.

## Matched 5 mm smoothing arm

```bash
bash code/functional_cluster_ica/start_smoothed.sh
```

This separate `nohup` submitter can be launched while the unsmoothed arm is
running. It records a waiting state, waits on the unsmoothed runner's lock,
then holds a shared lock to keep that source stable. The existing sm-0 runner,
its code and its output roots are not changed or restarted. No shared helper
implementation used by its frozen provenance is edited by this addition.

The matched arm uses 64 RAM-gated FEAT/extraction/task slots, six fixed ICA
slots and serialized automatic PCA. Check progress without large log dumps:

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/functional_cluster_ica/status_smoothed.py --brief
```

Every candidate sm-5 fit must match the unsmoothed arm's exact BOLD identity,
run-specific TR, complete design/contrasts, fitted original nuisance matrix,
geometry and other FEAT settings. Only smoothing/output destinations may
change; equivalent task-EV content is compared by hash. Candidates come from
saved original/accepted fit paths and the explicitly corrected TR control,
not a filesystem-wide search. Candidates that mismatch, including historical
UGR fits made from a different BOLD or the old TR, are not reused.

Only missing/mismatched fits are rebuilt in a new
`derivatives/insula_ica/functional-sm5-*` root, by cloning the active sm-0 FSF
and setting smoothing to 5 mm. A `feat_model` preflight requires the complete
temporal model to match before FEAT runs. FEAT smooths the **whole-brain BOLD
before ROI masking**, rather than blurring an already masked seed or smoothing
native whitened residuals after fitting. The FSL version must match the source
arm. No fMRIPrep or TEDANA is run. The two retained residual representations
then use the same algebra/extraction and physical masking as the sm-0 arm.

The ICA mask must remain **identical** to the unsmoothed arm. The code records
sm-5 validity and stops if any of those fixed voxels is unavailable; it never
silently erodes the mask or changes participants/runs to improve a match.
Balancing, input multiplier, concat, native variance normalization, 3/5/7/0
settings, deterministic starts, TR handling and convergence rules are shared.
Raw task-preserved ROI values and the functional/anatomical means are reused
from sm-0, with no new BOLD extraction or mean fitting. The new discovery maps
are spatially regressed onto that same raw task input before the established
pooled-UGR/full-task+nuisance standardized OLS. Thus the task comparison isolates
the smoothing used to discover components, keeping task measurement constant.

New records are `linux2-functional-sm5-*`. They publish selected/rejected
candidate reasons, FEAT output identities, coverage, ICA commands/maps and
five-bar figures. `smoothing_spatial_correspondence.json` matches sm-0/sm-5 maps
by absolute spatial correlation without task-based selection, retaining all
pairwise correlations and unmatched components. When corresponding sm-0 task
estimates are available, `smoothing_paired_task_comparisons.tsv` reports
sign-aligned, paired participant changes with descriptive pointwise CIs.
Component indices are not assumed to match. Original estimates are untouched.
Checkpoints reuse completed sm-5 FEATs, residual arrays, sessions, ICAs and
per-basis temporal fits on resume; failures and final/partial results publish.

Six additional tests verify that the generated FSF changes only allowed
settings; candidate reuse rejects different BOLD/TR/filter/design/confounds;
a completed FEAT checkpoint prevents any FSL refit; a mismatched regenerated
temporal model stops before FEAT starts; raw task caches require no
new image reads; and spatial permutation/sign matching precedes paired task
differences. Real Linux2 FEAT/convergence/QC checks remain execution-time checks.

## Whole-brain dual regression: fixed 3- and 5-component solutions

```bash
cd /ZPOOL/data/projects/rf1-betrayal
git pull --ff-only origin insula-ica
bash code/functional_cluster_ica/start_wholebrain.sh
```

The `nohup` launcher reuses the completed functional sm-0 and matched sm-5
solutions, both BOLD-scale and native-FILM residual representations. Eight
fixed bases contain 32 components in total. It has 32 RAM-gated run workers
(4 GiB reservation each plus shared 8 GiB headroom), then 32 single-threaded
`randomise` slots (2 GiB reservation each). The memory gate uses available
memory, so actual concurrency may be lower. There is no new FEAT, FILM, ICA,
fMRIPrep or TEDANA execution in this stage.

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/functional_cluster_ica/status_wholebrain.py --brief
```

Stage 1 regresses all component maps in a solution, plus a spatial intercept,
onto its **117 physically masked residual voxels**. The resulting timecourses
are centered and normalized to unit sample SD. Stage 2 regresses all these
timecourses together onto matching whole-brain residuals. BOLD-scale data use
the exact saved joint FILM fit: centered filtered BOLD minus the full design
times existing parameter estimates. The implementation applies this projection
in chunks, without storing another whole-brain 4D residual or refitting the
task/nuisance model. The native sensitivity uses the retained voxelwise
prewhitened `res4d`; its maps describe innovation-scale partial association,
which differs from unwhitened temporal covariance. Both dimensionalities share
each whole-brain read. Source hashes and retained output identities are checked.
Pilots include Trust, the pooled-pmod sub-11031 UGR run, and the corrected-TR
sub-10817 UGR run before expansion. Existing voxelwise full-model residuals
are retained; pooling pmods was an estimability repair for the separate raw
timecourse task model, not a new voxelwise fit here.

Runs are averaged within participant/task; the two tasks receive equal weight
within participant. The group input is **132 independent participant maps**.
Per-task group means are descriptive. Inference asks about pooled background
coupling, not task differences or social modulation of connectivity. The
whole-brain mask requires finite positive fitted variance and FEAT support in
all 980 run-arm entries. The entire original 127-voxel defining ROI is excluded
from group inference. No effect-based coverage selection is used.

Each of the 32 components gets 5,000 sign-flip permutations with 3D TFCE and
explicit positive and negative contrasts. The command is:

```text
randomise -i participant_maps.nii.gz -o dr -m inference_mask.nii.gz \
  -d group_design.mat -t group_design.con -n 5000 -T --seed=20261101
```

The design is a column of ones; contrasts are +1 and -1. Data are not demeaned
across participants (`-D` would remove the tested mean). FSL outputs `1-p`.
Reports distinguish per-map/tail TFCE FWE .05 from Bonferroni correction over
the two tails and all components within each basis (6 or 10 tests), and over
the entire requested set (64 tests). Within-basis and all-model signed maps
and figures are saved; all per-map/tail corrected probability maps remain
available. This is conservative across correlated solutions.

Anatomical component views and paired social task plots publish first. Then
group means publish before permutation inference. Final records include
TFCE maps, command arrays, mask coverage, anatomical target coefficients and
descriptive correspondence with the independent Smith et al. (2009) RSN10
reference. These latter summaries exclude the original ROI plus a predeclared
10-mm guard band, to reduce local spatial overlap; this extra guard band does
not alter the whole-brain inference mask. All template correlations are
reported, without forced network assignment. ICA signs and indices remain
local to each solution; do not equate IC numbers between arms.

Production maps stay under the separate ignored
`derivatives/insula_ica/functional-wholebrain-*` root. Small group/inference
maps, masks, plots, commands, tables and failures publish under
`code/functional_cluster_ica/records/linux2-functional-wholebrain-*` on
`insula-ica`. Repeating the same launcher verifies and reuses completed run,
participant, merge and permutation checkpoints. Partial attempts are preserved
separately. A lock prevents duplicate network runners and shared source locks
protect the completed functional arms. Failed publication is reported by the
status tool and retried on resume.

Five new tests verify whole-brain projection against direct residual regression
at 3/5 components and both residual scales, paired task contrasts, balanced
participant aggregation, frozen group designs and reuse without new external
permutations. A native FSL synthetic check validates the two signed one-sample
TFCE outputs and completed-checkpoint reuse. Its 100-permutation fixture is
not participant inference; production always uses at least 5,000 permutations.

Method references: [FSL dual regression](https://fsl.fmrib.ox.ac.uk/fsl/docs/resting_state/dualregression.html),
[FSL randomise](https://fsl.fmrib.ox.ac.uk/fsl/docs/statistics/randomise.html).

## Social versus nonsocial task characterization

Original five-bar estimates remain available. Additional four-condition views
show UGR nonsocial/social fairness slopes and Trust computer/stranger
reciprocate-minus-defect contrasts. The two primary paired differences are
**UGR social minus nonsocial** and **Trust stranger minus computer**; stranger
matches UGR's social partner more closely than friend. Friend-minus-computer,
friend-minus-stranger, and friend/stranger-average-minus-computer remain
secondary summaries in tables. Means and confidence intervals are computed
from each participant's contrast, preserving the paired design.

Do not average UGR slopes and Trust binary contrasts into a single coefficient,
or infer social specificity from one bar crossing zero and another not doing
so. Original coefficient directions remain: a negative social difference can
indicate stronger responses to social norm violations. Pointwise 95% CIs are
descriptive, without multiplicity-adjusted task inference; the functional ROI
was task selected. Testing whether connectivity itself is more social would
require a separate condition-interaction/PPI model.
