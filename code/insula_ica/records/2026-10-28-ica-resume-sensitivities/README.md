# ICA convergence recovery and separate sensitivities

The user authorized finishing the remaining ICA fits without redoing completed
work, examining voxelwise variance normalization, adding a separate unsmoothed
FEAT stream with 64 workers, and testing confounds-only/task-preserved ICA.
All streams retain N=132 and the accepted 241 Trust/249 UGR run identities.

## Linux2 commands

```bash
(
set -e
cd /ZPOOL/data/projects/rf1-betrayal
git pull --ff-only --quiet origin insula-ica
bash code/insula_ica/start_paired_group_ica_resume.sh
bash code/insula_ica/start_unsmoothed_feat.sh
bash code/insula_ica/start_task_preserved_ica.sh
)
```

Each starts a separate background runner and prints its PID/log. The launchers
have distinct locks; memory checks limit actual concurrency. Record publication
uses a shared file lock to serialize Git commits/pushes from concurrent streams.
Only each runner's specific new record is staged; unrelated work is preserved.

## What failed and the exact MELODIC setup

The original completed run record is
`linux2-paired-group-ica-20261028T215022Z-679a3e` (commit `311ae4d7`).
Underlying logs were published in
`linux2-paired-ica-diagnostics-20261028T233808Z-cf990d` (commit `a7c270f7`).

PCA and whitening completed for every failed fit. Symmetric FastICA then reached
2,000 iterations without convergence. Failures: joint-BOLD left dim 7;
joint-BOLD right dims 5/7; native-FILM right dims 5/7. The five trajectories have
large late changes, not just a residual barely above the convergence threshold.
There is no logged evidence of a memory kill, mask failure, or missing input.
Thirteen fixed fits passed, including all six bilateral fits. Six automatic fits
were not started because the original driver stopped at the fixed-fit gate.

Original command pattern (actual full paths are in the recorded plan/logs):

```bash
melodic -i sessions/RESIDUAL_MODE/ROI/inputs.txt -o ica/RESIDUAL_MODE/ROI/dim-DIM \
  --tr=1.615 --nobet -a concat -m masks/ROI_crop.nii.gz \
  --bgimage=reference_crop.nii.gz --report --Oall -d DIM \
  --disableMigp --vn --update_mask --seed=20261028 --maxit=2000 --maxrestart=10
```

This is temporal-concatenation group ICA over 264 participant/task blocks,
using classic group PCA. `--vn` switches off MELODIC's additional normalization;
the inputs already have within-run voxelwise unit-SD scaling and run weighting.
`--update_mask` switches off automatic mask updating. Separate physical 4D masking
is retained. In the inspected MELODIC source, the concat optimizer calls symmetric
FastICA; the restart option is used in the deflation path and did not provide ten
concat restarts. Thus explicit external starts are needed.

The resume changes only the iteration budget and random initialization: 10,000
iterations, seeds 20261028 through 20261032 in that order, the same `concat`,
`pow3` nonlinearity and 0.00005 convergence tolerance. Select the first converged,
validated fit; save every attempt. No tolerance relaxation, new objective, task
effect selection, dimensionality reduction, or component dropping is introduced.
This may resolve optimizer instability; it does not guarantee convergence.

Existing 13 fits are verified by their command, mask, input collection and map
identities, then adopted. The 490 source arrays and 1,584 masked inputs are checked,
not rebuilt. New attempts go under the existing root's `ica_retries/`; old failed
folders remain intact. `selected_ica_fits.json` identifies the actual selected
paths. The six automatic fits run independently of fixed retry success. Eleven
queue slots use 2 GiB reservations for fixed fits, 16 GiB for automatic fits, plus
8 GiB headroom. No fMRIPrep, FEAT or FILM is run by this resume. All 24 successful
fits are required before the final review/comparison phase.

## Mean, voxel count and native normalization

The fixed mask bypasses intensity-based mask creation, percent-BOLD conversion
is off, and MELODIC removes each voxel's temporal mean before PCA. A synthetic
264-input mean-zero versus exactly representable +64 control produced bit-identical
component maps locally. A +1024 control gave spatial correlations above
0.999999999998 with small amplitude differences from finite-precision input
rounding; adding large baselines needlessly reduces floating-point precision.
The production inputs stay mean-zero.

The shared masks retain all 391 left and 427 right voxels (818 bilateral).
Largest/smallest retained PCA eigenvalue ratios in the five failed fits are about
1.9–2.3: this does not indicate a numerical rank shortage. Spatial correlation
from smoothing can reduce effective spatial information. Bilateral convergence
is consistent with, but does not prove, a spatial sample-size issue. ICA component
identifiability and the number of biological subnetworks are not established by
voxel counts or convergence alone.

MELODIC's native normalization estimates variance after accounting for large
structured signals; it is not identical to the prior unit-SD scaling. The classic
joint-normalization source also has an absolute SD floor of 0.01 that can zero
channels. Our balanced inputs have small amplitudes, so simply removing `--vn`
can interact with that floor.

After the baseline attempts, run a separate paired pilot: right ROI, dim 5, both
residual types, with native MELODIC normalization enabled. Derived copies multiply
every value by the same factor 1000, preserving all relative run/voxel weights,
means and masks while avoiding the low-amplitude floor. Nothing is added to the
voxel mean. The native-normalization pilot has its own directory under
`native_vn_pilot/`, separate logs and comparisons. It never replaces a baseline
fit or changes only one member of a baseline pair. It is a numerical/method
sensitivity pilot, not a full 24-fit native-normalization series.

## Unsmoothed FEAT stream

`run_unsmoothed_feat.py` reads the active saved L1 FSF for each accepted run and
uses its exact BOLD path, task EVs and preserved complete nuisance matrix. It
copies that setup, changes `fmri(smooth)` from 5 to 0 and redirects output into
`derivatives/insula_ica/sm-0-feat-*`. Original models/contrasts, filtering,
prewhitening, timing and input identities remain required. Missing exact files
are recorded; available runs proceed without substituting smoothed filtered data.
No upstream preprocessing is launched. Nuisance and complete design identity,
residual geometry, finite/nonconstant coverage and retained intermediates are
checked. New task COPEs are expected to differ with smoothing and are not tested
for equality to sm-5 COPEs.

There are 64 FEAT worker queue slots with 1 GiB reservations plus 8 GiB headroom.
Call `feat` directly on the separate generated FSFs. Do not invoke original
wrappers: they may remove existing output folders and delete required residuals,
autocorrelations and (Trust) filtered data. All new intermediates are retained.
Use `--resume-record RECORD --execute --publish` to reuse an interrupted sm-0
stream. Completion still requires coverage/registration review before sm-0 ICA.

[Leech 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC6623935/) reports 6 mm
smoothing; the directly relevant restricted-ICA
[Leech 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621313/) used 5 mm and
explicitly included smoothing in its functional-connectivity preprocessing.
Thus the current 5 mm setting follows that precedent; sm-0 is a separate
sensitivity. Leech 2012 derived the restricted components from resting-state data
and characterized task responses afterward; task events did not drive that ICA.

## Task-preserved sensitivity

Use the existing sm-5 filtered BOLD to hold smoothing fixed. Build a matched OLS
pair using the saved complete FEAT design: (1) project out only its appended
nuisance block (the exact saved TEDANA/confound columns), and (2) project out the
entire task+nuisance design. Both start with voxelwise centered Y; no new temporal
filter or smoothing is added. The nuisance-only arm retains task variance outside
the nuisance span. Task-correlated nuisance variance is necessarily removed.
Do not subtract nuisance predictions from the old joint GLS fit as a substitute
for a genuine nuisance-only model.

Using matched OLS isolates task-regressor removal within this sensitivity and
avoids interpreting a GLS-to-OLS change as a task-removal effect. The original
GLS/innovation comparison remains unchanged. This is an explicit new sensitivity
model authorized by the user, not an additional cleanup of the baseline residuals.

Run only five components initially for each left/right/bilateral mask in each
matched arm: six fits. Use native MELODIC normalization consistently in both arms,
the same within-run SD/run weighting and global factor 1000, and physically masked
4D inputs. Per-run records quantify conditional task variance removed. No FEAT,
FILM or fMRIPrep is needed. New inputs/results go under a separate
`task-preserved-ols-*` root, with six component figures and descriptive map
comparisons. This is not an uncontrolled factorial expansion over smoothing,
residual representation, normalization and dimensionality.

See `validation.json` for local tests and source references. Linux2 outcomes,
including any remaining nonconvergences or missing sm-0 inputs, must be reviewed.
