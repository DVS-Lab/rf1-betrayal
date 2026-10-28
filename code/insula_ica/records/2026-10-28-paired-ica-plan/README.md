# Paired residual-scale insula ICA

The user accepted all three reviewed Trust flags, requested both residual
representations, required physically masked 4D inputs for each insula ROI, and
requested automatic dimensionality estimation in addition to the fixed dimensions.
No participants or accepted runs are removed: N=132, Trust 241 runs, UGR 249 runs.

## Linux2 launch

```bash
(
set -e
cd /ZPOOL/data/projects/rf1-betrayal
git pull --ff-only --quiet origin insula-ica
bash code/insula_ica/start_paired_group_ica.sh
)
```

The background launcher prints its PID and a short log command. Preparation uses
64 queue slots with RAM reservations (1 GiB for retained native residuals, 4 GiB
when a FILM replay is needed, plus 8 GiB headroom). The 18 fixed-dimensionality
fits run up to 18 at once, four threads each. The six automatic-dimensionality
fits follow in a separate queue: up to six slots, each requiring 16 GiB available
reservation plus shared 8 GiB headroom. The actual concurrency therefore depends
on available memory. No fMRIPrep, TEDANA, Warpkit or FEAT preprocessing is run.

Outputs stay under a fresh `derivatives/insula_ica/paired-group-ica-*` root.
Only its new `code/insula_ica/records/linux2-paired-group-ica-*` record is committed
and pushed, including failures. Original FEAT trees and unrelated working files
are preserved. Successful results remain `paired_ICA_complete_review_pending`;
this command does not run network connectivity or task-effect inference.

Use `python code/insula_ica/status_group_ica.py --brief` with the Linux2 Python
environment for a compact live report. It shows prepared inputs, validated fit
completion markers (fixed 18 plus automatic 6), active processes, automatic component
counts, and bounded error lines. Partial maps do not count as completed fits. A full
status snapshot is saved under ignored `code/insula_ica/output/paired-ica-status-*`.

## Accepted sources

- Trust base: `linux2-trust-regeneration-20261027T162428Z-ef53c5`.
- Replace 108 Trust run-2 fits with completed stock-filter outputs in
  `linux2-trust-stock-filter-20261028T203335Z-c78e1b-resume`; require same-run SBRef
  evidence. Retain reviewed sub-10720, sub-11326 and sub-11433 run 2, per the
  user's explicit instruction.
- Override both sub-10834 Trust runs with the accepted no-SBRef repair selected
  by `linux2-registration-fix-sub-10834-20261028T151154Z-acdf91`.
- UGR: the 249 included rows of the frozen Linux2 audit
  `output/linux2-audit-20261026T213907Z-dc0e41/audit/manifest.tsv`, retaining the
  audited filtered data and complete joint model. If a native residual is missing,
  replay only the exact logged FILM statistics command into a new directory.
  Require agreement with the retained model-defined COPE/ZSTATs and PEs.

Before computation, all runs must retain their full historical design/contrasts,
exact original nuisance matrix (including the saved TEDANA columns), timing,
compatible MNI grid, and coefficient images. Regenerated Trust maps need not
match the historical maps; that difference was previously accepted for this new
analysis. UGR sources must still match the original frozen inventory.

## Two residual representations

For each voxel, let Y be the retained `filtered_func_data`, X the complete saved
`design.mat`, and beta the existing joint FILM parameter estimates.

- `joint_bold`: `(Y - temporal_mean(Y)) - X beta`, on the original filtered-BOLD
  scale. There is no new fit, nuisance regression, temporal filter, or OLS cleanup.
  GLS residuals need not have zero ordinary Euclidean projection onto X.
- `native_film`: the native `stats/res4d.nii.gz`, on FILM's voxelwise prewhitened
  innovation scale. Its temporal covariance differs from original-scale residuals.

Five prespecified real-data pilots must verify native FILM residual reconstruction
and orthogonality using the corresponding voxelwise whitened design: an ordinary
Trust run, repaired sub-10834 run 1, retained sub-10720 run 2, and first/last UGR
runs. These checks precede the remaining preparation and ICA. Per-run diagnostics
also retain sampled lag-1 autocorrelations and between-representation temporal
correlations; these do not establish equivalence of the representations.

## Physical 4D masking, balancing, and dimensions

Use the existing Harvard–Oxford cortical maxprob-thr25 insula definition, split
by MNI x into left/right. Intersect fitted-mask, finite/nonconstant support across
all 490 runs and both representations. The common support is fixed between scales
and dimensions. Record per-run and common coverage; stop below the existing 80%
common-coverage gate without automatically dropping runs.

Each participant/task contributes one PCA input block (264 blocks per fit).
Keep every accepted run/frame. Demean and standardize each voxel within each run,
then weight that run by `1/sqrt(number_of_runs * run_frames)` before concatenating
within participant/task. Each run has equal squared energy within its session.
Disable MELODIC's additional variance normalization and automatic mask updating.
No extra temporal filter is added.

**Write separate left, right, and bilateral 4D NIfTIs for each representation.**
Every voxel outside the specific ROI is exactly zero. Read and verify the saved
4D inputs; save hashes and validation results, and recheck identities before
each fit. Supply these physically masked files and the corresponding explicit
mask to MELODIC. A crop of the common bilateral bounding box preserves MNI voxel
coordinates without interpolation. Component maps are also padded back onto the
full original MNI grid.

For both representations and each ROI, run **5 primary**, **3/7 sensitivity**,
and **0 (automatic estimation)**: 24 fits. Use classic temporal-concatenation
group PCA (`--disableMigp`), fixed seed 20261028, and retain all components.
Record the actual selected component count for every automatic fit. The counts
may differ between representations. Missing maps/nonconvergence fails the fit
and is recorded; no partial fit is treated as completed.

Compare the two scales within each ROI/dimension using all-pair spatial
correlations, one-to-one maximum absolute-correlation matching, and centered
subspace angles. Allow rectangular matching for automatic fits, retaining unmatched
component identities. Compare all dimensionalities against the 5-component basis.
Signs and component indices are arbitrary. Do not choose dimensionality or retain
components based on downstream task effects. Save all component figures, paginated
for large automatic solutions, and complete numerical comparisons for review.

## Validation and resuming

See `local-validation.json`. Unit tests include actual recorded Trust source
selection, the 490-run orchestration with mocked fits, physical saved-volume
mask checks, identity/tamper checks, rectangular matching, and recovery boundaries.
Local real FSL tests verify the FILM residual definitions and fixed/automatic
MELODIC operation on synthetic data. Development FSL is 6.0.7.19; Linux2's actual
version and real-data pilot outcomes are recorded at runtime. Synthetic success
does not imply that the Linux2 real-data pilots or production fits have passed.

The launcher prevents duplicate runs. For a failed run, inspect its record first.
Validated source arrays, completed UGR recoveries, balanced sessions and completed
ICA fits can be reused through `run_group_ica.py --resume-record RECORD --publish`.
Partial files without completion markers are preserved and require inspection;
the script never silently overwrites them. Resume validates the same sources,
masks, settings and saved outputs. Git publication stages only the specific new
record and preserves unrelated local/staged work.
