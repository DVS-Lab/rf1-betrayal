# Insula ICA from preserved first-level FEAT residuals

This workflow follows the revised 2026-10-26 brief. **The earlier workflow that
added a new nuisance-cleaning model before ICA is retired.** Old configs are
rejected. Existing first-level FEAT task+nuisance fits supply the residual data;
original task COPEs supply a separate task-preserved analysis branch.

The earlier user addition remains in force: **left, right and bilateral insula**,
each with one shared TG+UG basis, five primary components and 3/7 spatial checks.
Retain every component. The 30 primary ROI × component × task tests form one
Bonferroni family. Component numbers do not match across ROIs/dimensionalities.
No new gPPI or behavioral search is included.

## Current evidence and limits

Jacob's corrected `master_subject_table_full (1).csv` has 147 rows and exactly
132 `include_final=TRUE` identities, matching `code/sublist_n132.txt`. For that
cohort, the saved full-run QC permits:

| Task | Run 1 | Run 2 | Total accepted runs |
|---|---:|---:|---:|
| Trust | 131 | 110 | 241 |
| UGR | 132 | 117 | 249 |

These are **QC-eligible counts, not verified linux2 image availability**. The
private CSV and ages are not committed. The startup helper can use a saved
linux2 copy or restore the supplied cohort-only snapshot into ignored output.
It rejects the superseded 137-included table. Keep the manuscript list unchanged.

Source inspection of `L1stats-trust.sh` and `L1stats-ugr.sh` confirms both delete
`stats/res4d.nii.gz`, `stats/corrections.nii.gz`, and `stats/threshac1.nii.gz`.
Trust also deletes `filtered_func_data.nii.gz`. Both wrappers can skip completed
models or remove an existing FEAT directory. **Never run those wrappers to
recover these inputs.** The new scripts do not call them or change them.

The model templates use TR=1.615 s, smoothing=5 mm, zero volume deletion,
disabled FEAT temporal high-pass, and **FILM prewhitening enabled**. Per-run
saved settings and designs, not the templates alone, are authoritative.
The old `fmriprep-24` images may be missing. Neither the old images nor historical
FEAT trees are available on the development Mac, so no real-data reconstruction
or historical comparison has yet been performed.

## Trust preprocessing rebuild after the linux2 inventory

The user's 2026-10-26 linux2 inventory found a FILM replay route for all 249
eligible UGR runs, but no residual, filtered data, or source BOLD at the recorded
paths for all 241 eligible Trust runs. The subsequent direct RF1 BIDS check found
four raw magnitude echoes plus JSONs for all 241 Trust runs and session-01 T1w
files for all 132 participants. No other project needs to be searched.

`rebuild_trust.py` uses the archived cluster configuration from
[`a19c46ef:code/fmriprep24-hpc.sh`](https://github.com/DVS-Lab/rf1-sra-linux2/blob/a19c46ef/code/fmriprep24-hpc.sh)
and its `ses-01.json`: fMRIPrep **24.1.1**, multi-echo outputs, T1w and
MNI152NLin6Asym spaces, no FreeSurfer reconstruction. It restricts functional
processing to accepted session-01 Trust runs and adds explicit resource caps.
The archived command is evidence for the recovery configuration; it does not
prove that today's raw BIDS data, metadata, fieldmaps, or template cache are
identical to the inputs used for every historical analysis.

**No standalone TEDANA rerun is needed for historical FEAT recovery when its
original combined confounds remain.** The inventory reached its missing-BOLD
checks only after confirming a surviving external or FEAT-preserved matrix.
The rebuild checks every run's matrix again, preserves its original bytes, and
uses it unchanged, including the original rejected TEDANA component time
courses. The new fMRIPrep confounds and any new ICA decomposition are not
substituted into that historical model. Re-running standalone TEDANA would
produce a new nuisance fit, rather than recover the old one.

Start the three-run, two-participant pilot (one participant with two accepted
runs and one with a single accepted run):

```bash
INSULA_AUDIT_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
INSULA_AUDIT_DIR=/ZPOOL/data/projects/rf1-betrayal/code/insula_ica/output/linux2-audit-20261026T213907Z-dc0e41/audit
"$INSULA_AUDIT_PYTHON" code/insula_ica/rebuild_trust.py \
  --audit-dir "$INSULA_AUDIT_DIR" --pilot --fit --fetch-image --execute --jobs 2
```

Remove `--execute` to inspect the plan first. `--fetch-image` pulls the official
`nipreps/fmriprep:24.1.1` container into the isolated rebuild root if the exact
`/ZPOOL/data/tools/fmriprep-24.1.1.simg` is absent; `--image /path/to/image` can
use an existing cluster container. The runtime version is checked and the
container SHA256 recorded. No existing environment or current fMRIPrep wrapper
is changed. Defaults are two concurrent participants for the pilot and **eight
for production**, each with 14 processes and a 28,000 MB memory cap. Explicit
`--jobs` overrides concurrency; `--nprocs` and `--mem-mb` set per-subject limits.
Eight concurrent subjects therefore have 112 process slots and 224,000 MB of
combined configured memory caps.

The runner preflights raw echo volume counts/TR against each saved FSF, preserves
the full combined nuisance matrix and original FSF, and records raw image
size/mtime and JSON hashes. New outputs go under:

```text
derivatives/insula_ica/trust-rebuild-24-<audit-session>/
  fmriprep-24/           # new BOLD; original data remain untouched
  original_confounds/   # byte-preserved historical combined matrices
  historical_fsf/
  filters/
  feat_reconstruction/  # with --fit; retained res4d/filtered data/threshac1
  sub-*.fmriprep.log
```

Work directories are isolated under `/ZPOOL/data/scratch/$USER/rf1-betrayal-trust24/`.
BIDS is bound read-only. Completed fMRIPrep participants can resume using their
completion marker; interrupted work remains available to the container. An
execution lock prevents duplicate launches. A partial FEAT is preserved and
requires inspection; the runner never invokes the original deleting wrappers.

The archived cluster and current linux2 Warpkit scripts both export the products
consumed by fMRIPrep into `bids/sub-*/ses-01/fmap/` as
`sub-*_ses-01_acq-trust_run-*_fieldmap.nii.gz`, `_magnitude.nii.gz`, and their
JSON sidecars. They do not overwrite the raw `func/` BOLD echoes. The runner
checks these exact four files, one-volume fieldmap/compatible magnitude grid,
`Units=Hz`, and `IntendedFor` or matched `B0FieldIdentifier`/`B0FieldSource`
associations to all four raw magnitude echoes. It records all 241 run results
in `warpkit_inventory.json`, selecting ready participants for the pilot and
requiring every selected run to be ready before execution. A `.warpkit_done`
flag and the large Warpkit intermediate files are not required when those
usable BIDS products remain. The exact retained 4D fieldmap location is also
reported if a BIDS export is missing; the runner never invokes Warpkit or
silently runs fMRIPrep without the required distortion-correction inputs.

With `--fit`, only three saved FSF locations change: input BOLD, the preserved
confound file, and output directory. All original task/nuisance settings remain.
The runner compares the complete design/contrast matrices and every historical
COPE/ZSTAT (atol=rtol=1e-4). A difference stops expansion; it is not accepted as
an equivalent historical reconstruction. The regenerated BOLD must match the
historical spatial grid, TR and time count. Full residual-scale/coverage QC is
still required before ICA; a successful fit comparison is not that final QC.

After the pilot passes, launch the full rebuild with eight concurrent subjects:

```bash
"$INSULA_AUDIT_PYTHON" code/insula_ica/rebuild_trust.py \
  --audit-dir "$INSULA_AUDIT_DIR" --fit --fetch-image --execute --jobs 8
```

Full execution requires the matching `pilot-verified.json`.
The runner does not automatically feed new Trust residuals into ICA; combining
their reconstruction manifest with the existing UGR replay route remains the
next step after real-data historical comparisons pass.

Run filters use JSON string values (`"run": ["1", "2"]`), because the 24.1.1
configuration deserializer crashes on integer run values. Before launching any
participants, the runner exercises the installed container's own JSON/filter
converter functions on every selected filter. Failed startup attempts can be
restarted with the same command/output root; corrected filters are regenerated,
the downloaded container and work directories are retained, and participant
logs append a timestamped attempt marker.

Container identity compares the image SHA256 and parsed application version,
not the entire `--version` stdout/stderr stream. Update notices, proxy warnings,
and transient network messages are appended to `container_observations.jsonl`.
Legacy records storing the complete command output are normalized only after
the SHA256 and parsed version match; the original bytes remain in
`container.legacy.json`. A genuinely different image or version still blocks,
with the original `container.json` left intact. Do not delete the identity record
or change the output root merely to bypass a notice-only difference.

Development checks: thirteen reconstruction safety tests plus ten existing numerical
tests passed. A 241-run synthetic inventory exercised the full raw preflight,
three-run pilot selection, byte-preserved confounds and isolated output plan.
The exact upstream 24.1.1 converter functions reproduced the integer-run error
and accepted the corrected filters in a source-level check with a Query enum
fixture. The actual container preflight runs on linux2.
No fMRIPrep container or real RF1 images were executed on the Mac.

## Important residual-scale gate

The installed FSL 6.0.7.19 source (`src/fsl-film/film_gls.cc`) prewhitens each
voxel's data/design, computes residuals in that transformed space, and writes
them as `stats/res4d`. An actual synthetic `film_gls` run confirmed this:
prewhitened-data minus the matching whitened design × PEs reproduced `res4d`
(relative error <1e-4), while original-scale `Y-X*beta` did not.

Therefore **a FEAT filename alone does not establish an unwhitened background
BOLD series**. Voxel-dependent whitening changes the temporal covariance being
characterized. The scripts verify the native product and stop before ICA unless
its scale is acceptable under `residual_scale_review`.

Default `accept_voxelwise_prewhitened=false` deliberately leaves that decision
open for the real pilot. Set it to true only after explicit scientific review,
with a written `reason`, if innovation-scale residual networks are the intended
estimand. Do not interpret those maps as original-scale BOLD connectivity.
No inverse whitening, alternate residual reconstruction, extra regression or
additional temporal filtering is silently applied.

If unwhitened full-model residual connectivity is required, the alternative to
review is the **original joint FILM decomposition**: centered
`filtered_func_data - design.mat * original_FILM_PEs`. This uses the same fitted
model without fitting another nuisance regression. It is proposed here, not
automatically substituted for native `res4d`. It requires retained/recreated
filtered data and a reviewed scale/orthogonality validation.

## First linux2 operation

Use the upstream-documented Python interpreter:

```bash
INSULA_AUDIT_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
"$INSULA_AUDIT_PYTHON" code/insula_ica/audit_linux2.py --qc /path/to/corrected/master_subject_table_full.csv
```

`audit_linux2.py` prints its actual Python and FSL versions, searches when no
`--qc` is supplied, checks the corrected cohort identity, creates fresh output
and config, and **runs inventory only**. `--qc-snapshot-stdin` reads the supplied
gzip/base64 cohort-only CSV; if no matching saved CSV exists, it restores that
snapshot locally. No shared Python environment is changed and no FEAT is run.

The printed paths have this structure (the timestamp/token is generated once):

```text
<betrayal>/code/insula_ica/output/linux2-audit-<timestamp-token>/
  config.json
  audit/report.md
  audit/manifest.tsv
  audit/<run-key>/reconstruction.fsf
<betrayal>/derivatives/insula_ica/feat-residuals-<timestamp-token>/
  feat_reconstruction/<run-key>.feat/
```

Both trees are ignored. The manifest gives the exact original FSF, historical
FEAT, BOLD and confound paths, whether BOLD/filtered data remain, actual model
settings, recovery route, target directory, and exclusions. Raw paths are read
from the generated FSF. No newer derivative is substituted.

You can also copy `config.example.json`, set `master_qc`, `output_root` and
`derivatives_root`, and call `insula_ica.py --config FILE audit`. A new run needs
fresh output/derivative roots. Relative config paths resolve against this repo.
Legacy pre-ICA-cleaning configs cannot be used.

## Recovery routes

The inventory chooses a route for each originally accepted run:

1. **Reuse existing `stats/res4d`** if present, followed by full residual QC. No
   requirement to restore raw BOLD solely to read a valid retained residual.
2. **Full FEAT reconstruction** if original BOLD, confounds and EVs remain. Copy
   the original generated FSF, change **only `fmri(outputdir)`**, and execute
   `feat <copied.fsf>` in the new reconstruction root. Preserve all intermediates.
3. **Stats-only FILM replay** if raw BOLD is gone but original FEAT
   `filtered_func_data`, full design/contrasts and the exact FILM command remain.
   Parse that command from the original log, retain fit options, use the same
   input/design, and change only the FILM results destination to new `stats/`.
   Copy small model/reference files into the new FEAT directory. Do not rerun
   smoothing, filtering, TEDANA or source preprocessing. A missing/ambiguous log,
   voxelwise EV model or unsupported command is a blocker, not an invitation to
   assemble an approximate command.
4. **Blocked** if no exact route exists. Report which original inputs are absent.
   If Trust lost both its source BOLD and filtered data, its residuals cannot be
   reproduced from COPEs alone. Restore historical inputs if available, or use
   the Trust preprocessing rebuild above and test its historical equivalence.

Reconstruction never writes in the historical FEAT directory. Completed new
residuals are not rerun. Partial new directories are preserved and require
inspection; there is no `rm -rf` recovery path or overwrite switch.

## Pilot then full reconstruction

Use the `config.json` path printed by the inventory:

```bash
INSULA_CONFIG=/exact/printed/path/config.json
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope pilot reconstruct
```

This is a **dry run**, writing the exact argv/cwd/output plan. It selects up to
five representative runs spanning tasks, recovery routes and model rank.
Missing runs elsewhere do not prevent a runnable pilot from being inspected.
After checking that every destination is in the new root:

```bash
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope pilot --execute reconstruct
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope pilot verify
```

`verify` compares **every historical COPE and ZSTAT**, full design and contrast
matrices, geometry/time count, finite residual coverage, and fitted-regressor
projection using the appropriate voxelwise whitened design. When filtered input
and PEs are present, it also reconstructs the native residual scale numerically.
The diagnostic FFT design transform is checked against real FILM; it is never
applied as another cleaning step to the analysis residuals. Do not compare
prewhitened residuals with unwhitened EVs and call the difference failed denoising.

Default voxelwise COPE/ZSTAT comparison tolerances are atol=rtol=1e-4. A material
difference stops the run; do not loosen tolerances to conceal a version/model
change. Original nuisance matrices and all model columns are retained. Complete
commands and source-file identities are saved. A surviving external confound
matrix must agree with FEAT's retained `confoundevs.txt` when both exist.

A successful pilot historical comparison is required for full reconstruction:

```bash
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope full reconstruct
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope full --execute reconstruct
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope full verify
```

Full reconstruction recovers the remaining inputs; ICA still requires successful
full QC and the scale review. FEAT may submit jobs on a configured scheduler;
if files are not ready when the command returns, inspect scheduler/log status
and later run `verify` rather than rerunning into that partial destination.

## ICA, residual connectivity and task-preserved branch

After all original accepted runs pass comparison/coverage checks and residual
scale is approved, `prepare` constructs fixed anatomical/valid-data masks and
references the native residual files **without a new cleaning/scaling operation**.
The HO cortical maxprob thr25 XML identifies Insular Cortex as label 2; left and
right are defined by MNI coordinates, with nearest-neighbor atlas resampling.
Intersect with common finite nonconstant residual coverage; save loss masks and
counts. Default retained ROI coverage gate: 80%; inspect anterior/ventral loss.

```bash
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" prepare
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope pilot ica
```

ICA pools all accepted TG+UG runs in one shared basis per ROI. One concatenated
input per participant/task balances the number of subject-level PCA blocks;
within that input, run/voxel variance scaling and `1/sqrt(n_runs*T_run)` weights
balance runs. All original residual frames remain; no new censoring is applied.
Standard MELODIC individual PCA/group concatenation is used. It performs ordinary
demeaning/spatial ICA scaling but **no additional temporal high-pass, band-pass,
new nuisance regression, or smoothing**. `--vn` disables its voxel variance
normalization, `--disableMigp` uses the explicit individual-PCA route, and
`--update_mask` disables automatic mask changes; output mask equality is checked.
Installed help and complete effective commands are saved. Freeze the stochastic
solutions rather than claim bitwise reproducibility from an unsupported seed.

Original maps are retained, with oriented copies whose largest absolute loading
is positive. No task coefficient or network label determines the sign. Save all
3/5/7 spatial correspondences; do not select a dimension based on task effects.

Set `ugr_pmods_increase_with_fairness` and `ugr_coding_evidence` after inspecting
stored pmod coding. Then run:

```bash
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope pilot regress
```

- **Connectivity stage 1:** one simultaneous spatial model with all five frozen
  maps and a spatial intercept into insula residuals. Save overlapping-component
  partial timecourses and the spatial projection operator.
- **Connectivity stage 2:** simultaneously regress those five TCs into the **same
  native whole-brain residual series**, after documenting TC SD scaling. Do not
  append task/confound EVs or use preprocessed BOLD as the target. Maps are
  descriptive partial/background coefficients, not bivariate seed correlations.
- **Task response:** apply the **same frozen multivariate spatial operator** to
  each run's **original task COPE**, separately from the residual timecourses.
  This is linear spatial projection of a task-preserved estimate, not recovery
  of modeled task signal from `res4d`. Negate Trust contrast 10 `rec-def`; use
  UGR contrast 11 across all partner/endowment pmods and reverse it only when
  increasing pmod amplitude means increasing fairness. Preserve original contrast
  sum/units. Task coefficients are not divided by residual TC SD.

No independent-voxel assumption is used to manufacture a spatially weighted
COPE variance. Within-participant projection SEs are left unestimated; equal-run
participant means enter group one-sample tests/CIs using empirical between-person
variation. All 30 primary tests receive Bonferroni correction. Task units differ;
no raw TG–UG interaction or matched cross-ROI component comparison is claimed.

Supply independent 3D MNI network templates (`name`, `path`, `citation`) before
`summarize`. Only grid interpolation is performed; templates must already be in
MNI space. Correlate outside the defining insula ROI, report weak/ambiguous
matches, and retain unlabeled/artifact components. Group connectivity averages
runs within participant/task, then tasks, then participants equally.

```bash
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope pilot summarize
# Review pilot maps, scales, source timing/sign and model/QC reports first.
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope full ica
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope full regress
"$INSULA_AUDIT_PYTHON" code/insula_ica/insula_ica.py --config "$INSULA_CONFIG" --scope full summarize
```

Outputs include all 15 component anatomies, descriptive network maps/template
correspondence, all 30 task estimates/N/SE/CI/corrected-p, two composite figures,
spatial stability, and `results_summary.md`. Potential component differences
remain exploratory; different maps alone are not a formal difference test.
The result is not independent validation or a unique betrayal/fairness encoding
claim, and does not test whether the original sample-defined aIns–vmPFC gPPI
survives an independent seed. Show null and ambiguous outcomes.

## Verification and provenance

```bash
"$INSULA_AUDIT_PYTHON" -m unittest discover -s code/insula_ica -p 'test_*.py' -v
"$INSULA_AUDIT_PYTHON" code/insula_ica/test_insula_ica.py --integration /absolute/new/scratch/path
```

The synthetic integration uses actual installed `feat_model`/`film_gls`/MELODIC,
creates historical fits, deletes their residuals, omits the raw BOLD, dry-runs and
replays stored FILM commands into separate roots, compares all original results,
checks native residual scale, and exercises all three ROIs and task/summary
branches. Historical files must remain byte-identical. A separate prewhitened
FILM probe validates the diagnostic whitening convention. No RF1 pilot result
is implied by these software tests.

This branch started at betrayal main `affb3b30ad9b911232a99e32df6459403c84d4ea`.
Original code/template/FEAT paths, input size/mtime, lightweight hashes, code
hashes, software and generated commands are frozen in each inventory. Existing
inputs cannot change between stages. Full private 4D hashes are not calculated.
Retain partial/failed outputs for diagnosis. Do not commit source imaging,
reconstructed 4D files, full private QC tables or ages.

Methodological sources:
[Leech 2012](https://doi.org/10.1523/JNEUROSCI.3689-11.2012),
[Leech 2011](https://doi.org/10.1523/JNEUROSCI.5626-10.2011),
[Lindquist 2019](https://doi.org/10.1002/hbm.24528), and the installed FILM source
(`film_gls.cc`, `AutoCorrEstimator.cc`, `glimGls.cc`). The supplied revised brief
controls this workflow. In particular, no additional sequential nuisance filters
are introduced after the joint FEAT task+nuisance fit.

Development checks on 2026-10-26: ten numerical/safety tests passed. The real-FSL
synthetic residual-route integration reconstructed six deleted run residuals from
retained filtered data with missing raw BOLD; all 210 historical COPE/ZSTAT maps
matched exactly and original files remained byte-identical. It completed all
3/5/7 decompositions for three ROIs in pilot/full scopes, residual connectivity,
original-COPE projection and the 30-test summaries/figures. The separate native
prewhitened FILM probe matched the whitened joint decomposition to 1.4e-5 relative
error. Full FEAT prestats reconstruction still requires the real linux2 pilot;
these synthetic results do not establish RF1 input availability or equivalence.
