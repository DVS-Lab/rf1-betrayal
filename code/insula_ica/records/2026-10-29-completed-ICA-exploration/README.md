# Explore completed components while automatic ICA continues

All 18 fixed-dimensionality unsmoothed fits completed according to the live
linux2 status. The source runner publishes at its end, so convergence alone did
not expose new anatomy/network/task results on GitHub. The user authorized review
and further analysis now while the remaining automatic fits run.

```bash
bash code/insula_ica/start_ica_exploration.sh
```

The exploration uses a read-only snapshot of completed fixed fits in the latest
corrected unsmoothed record. All six primary five-component solutions (both
residual representations × left/right/bilateral) must be complete. Completion
map identities, commands, fixed masks and normalization policy are checked. No
fit is selected by task/connectivity findings. No ICA, FILM, FEAT or preprocessing
job is started or stopped. The automatic ICA submitter remains independent.

## Published checkpoints

1. **Map review:** all completed fixed-fit component figures, six primary group
   map NIfTIs, peak/absolute-loading-center coordinates and bilateral energy
   asymmetry, 3/5/7 spatial comparisons, and paired residual-scale comparisons.
   These are signed loadings, not disjoint anatomical parcels or network labels.
2. **Task review:** all 60 task projections, participant means, two-sided tests,
   pointwise CIs and figures, published before whole-brain connectivity finishes.
3. **Connectivity review:** group maps and component/connectivity figures, an
   explicit common coverage mask, anatomical target means and signed rankings.
   These maps are descriptive background coefficients, not voxelwise group
   inference, task-modulated connectivity, or independently assigned networks.

## Connectivity calculation

The primary spatial model contains all five centered component maps and a spatial
intercept inside the defining ROI. Its pseudoinverse supplies the same projection
operator for residual timecourses and task COPEs. Timecourses are scaled to unit
sample SD and all five enter the temporal model together. A conditioning/rank
check rejects an unestimable basis.

Both stages use the same residual representation. Native FILM uses the retained
whole-brain `res4d`, so it describes innovation-scale coupling. Joint-BOLD uses
the existing full joint fit on the original filtered-BOLD scale. For each spatial
chunk the calculation is

```
A @ (Y - mean_time(Y)).T - (A @ X) @ B.T
```

where `A` is the pseudoinverse of standardized component timecourses, `X` is the
saved complete design and `B` the saved voxelwise GLS coefficients. This equals
projection of the same full-model GLS residual used to derive the ICA inputs;
it does not fit another nuisance model or apply another temporal filter. The
insula arrays, source data/model identities and corrected task COPE identities
are checked. No task/nuisance EVs are appended again to these already cleaned
residuals. This implements the two multiple-regression stages described in
[FSL's dual-regression documentation](https://fsl.fmrib.ox.ac.uk/fsl/docs/resting_state/dualregression.html),
using the previously agreed joint-residual inputs rather than sequential denoising.

Coverage is the intersection of fitted masks and finite positive sigma over all
490 runs. Missing coverage is distinct from zero connectivity. The defining seed
ROI is excluded from connectivity displays and anatomical target summaries.
Contralateral insula can remain a target for unilateral seeds. Cortical atlas
names describe regional targets; no salience/default-mode/network label is
automatically assigned. Signed coefficient magnitudes across the two residual
scales are not treated as interchangeable.

Runs are averaged equally within participant/task. Group maps average participants
equally (N=132); pooled maps also give each task equal weight. Run/participant
maps and timecourses are retained under a separate ignored derivatives root.
Only group spatial maps, coverage, review figures, derived tables and provenance
are copied into the Git record. Display cutoffs are descriptive, not significance
thresholds. Formal voxelwise inference and reference-network matching remain
separate analyses.

## Task calculation and interpretation

Already task-cleaned residuals cannot supply the original evoked effect simply
by regressing task EVs into them again. Instead the same frozen spatial weights
project each run's current sm-0 task COPE, including the corrected 1.710-second
UGR fit. Directions are kept **as stored**: Trust `rec-def`, UGR offer pmod.
No unverified increasing-unfairness sign flip is applied. The relevant EV coding
should be checked before attaching a directional psychological interpretation.

Runs are averaged within participant/task before a one-sample group test; runs
are not independent subjects. All 60 tests (two residual bases × three ROIs × five
components × two tasks) are reported, with Bonferroni-60 as the combined family.
Within-scale Bonferroni-30 values are also retained. Component indices do not
identify matched components across ROIs or scales. Maps were estimated in this
cohort, so all tests are exploratory, not independent confirmation, unique
psychological encoding, a formal between-task/hemisphere difference, or validation
of the original PPI. Null and ambiguous findings remain visible.

## Execution and validation

Eight whole-brain readers use four-GiB RAM reservations plus eight-GiB headroom
while automatic PCA continues. No active source output is modified. Completed
exploration runs have frozen signatures/output hashes for explicit resume via
`--resume-record`; unknown partial outputs are preserved for inspection. At least
32 GiB free disk is required for retained maps. Git publication is serialized
with other analysis publishers. Publication failures are recorded locally.

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/insula_ica/status_ica_exploration.py --brief
```

Seven local tests passed, covering exact GLS-residual projection, multivariate
spatial/temporal operators, temporal SD normalization, run-before-participant
aggregation, all 60 corrected tests, rejection of incomplete cohorts/primary
bases, read-only snapshots while source ICA is running, and real NIfTI whole-brain
outputs for all six bases matching direct residual regression. A real NIfTI
figure/export smoke check produced six primary figures/maps and 30 anatomy rows
at `/private/tmp/rf1-ICA-exploration-review-8ujgt9mh`. Shell syntax checks passed.
These are software checks, not findings about RF1 connectivity/task effects.
