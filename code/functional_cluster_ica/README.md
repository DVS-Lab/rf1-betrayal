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
