# Run-specific TR correction and project integration

The user authorized a 1.710-second reanalysis of
`sub-10817_ses-01_task-ugr_run-1`, retaining investigation of historical source
provenance. All four current raw JSONs and headers agree at 1.710 seconds; the
historical FEAT model used 1.615 seconds. The existing sm-0 recovery completed
489/490 FEATs. No raw metadata/header edit is authorized or needed.

## Linux2 correction

```bash
bash code/insula_ica/start_ugr_timing_repair.sh
```

The job audits raw echo TR/header consistency and original modeled TR across all
490 accepted runs, including the participant's two Trust runs. Selected protocol
fields are retained to help distinguish a run/sequence variation from a
participant-wide issue. No assumption that all of one participant's tasks share
TR is made. Other discrepancies are recorded rather than automatically repaired
or excluded.

Only the authorized UGR case enters fMRIPrep 24.1.1, with its verified cached
anatomy, existing Warpkit fieldmaps, unchanged stock filter and a one-run BIDS
view. Current raw image identities and JSON hashes must match the reviewed timing
evidence. The generated BOLD must retain the measured 1.710-second TR.

The existing seconds-based task event files are checked for finite valid values
and fit within the corrected acquisition duration. Their onsets/durations are
not multiplied by a TR ratio. This does not independently establish the scanner
zero offset or historical DICOM provenance; those remain investigation tasks.
The original saved TEDANA/nuisance matrix is retained exactly.

Three FEAT fits follow, in separate locations:

| Fit | BOLD | Smoothing | Model TR | Purpose |
| --- | --- | ---: | ---: | --- |
| Corrected sm-0 | New 24.1.1, TR 1.710 | 0 mm | 1.710 | Complete the missing unsmoothed input |
| Wrong-TR control | Same new BOLD | 5 mm | 1.615 | Diagnostic only |
| Corrected sm-5 | Same new BOLD | 5 mm | 1.710 | Corrected analysis candidate and matched comparison |

`feat_model` regenerates the task design at each TR. Corrected fits must have a
changed design; changing only an output header is explicitly rejected. Contrast
definitions and nuisance values must stay unchanged. The sm-5 wrong-TR control
must match the original temporal design. The same-BOLD comparison therefore
isolates the TR/model-sampling effect from new preprocessing. All COPE/ZSTAT
contrasts, including UGR contrast 11, are compared in the shared brain mask.

The corrected sm-0 output fills the previously absent canonical target in the
separate unsmoothed analysis root. Existing 489 successes are validated and reused.
Original sm-5 manuscript fits and old ICA inputs are preserved. The combined
sm-0 record explicitly labels the regenerated design and timing correction.
Outputs, event/source identities, Trust checks, cohort audit, matched statistics
and any failures publish to `insula-ica`. `corrected-analysis-input.json` identifies
the sm-5 replacement candidate; it is not silently substituted into old analyses.

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/insula_ica/status_ugr_timing_repair.py --brief
```

## Project-level fix prepared on this branch

`code/run_timing.py` validates all four raw JSON/header TRs and volume counts,
including NIfTI time units, and checks pre-FEAT BOLD timing against them. It rejects
inconsistent inputs rather than letting FEAT mask a discrepancy by relabeling its
output header. Its `stamp` command replaces only the generated FSF's numeric TR
and records a timing sidecar. Task-event files remain in seconds.

Both `L1stats-ugr.sh` and `L1stats-trust.sh` call this validation before their
existing-output skip. They refuse to silently reuse a fit with a different TR and
stamp the validated per-run TR before either activation or PPI FEAT execution.
They use `$FSLDIR/bin/fslpython`, overridable with `RF1_TIMING_PYTHON`. This is
normal FSL Python with nibabel/numpy, not an assumption about the user's base
conda environment. Upstream fMRIPrep scripts, raw metadata, task event generators,
and stored historical designs are not modified.

This prepares a project integration change for review/merge. It does **not**
retroactively repair existing activation/PPI/group results. After the timing audit
and matched case comparison, the project repair should:

1. Resolve original acquisition provenance and task scanner-zero timing; inspect
   any additional cohort discrepancies identified by the audit.
2. Promote versioned corrected first-level activation fits with explicit input
   paths. Recompute affected PPI/deconvolution regressors using the correct run TR;
   do not reuse stale physiological/task design caches.
3. Recompute dependent participant L2 and group L3 activation/PPI results, summary
   extractions and RSA inputs as applicable, recording before/after inference.
4. Replace the affected insula residual input through an explicit corrected input
   registry and regenerate dependent ICA preparations/fits. Existing paired/group
   ICA code still assumes 1.615-second session headers and must be updated to
   preserve mixed acquisition TRs before using the new registry. Do not relabel
   the 1.710-second data, interpolate it silently, or report frequency axes as if
   every block shared one TR.

Successful local FEAT correction is not a claim that the paper's group statistics
have been repaired. Those downstream reruns remain tracked work.

## Validation

38 targeted tests passed, including acquisition/header consensus and time-unit
conversion, both project launcher paths, corrected design acceptance, rejection
of a header-only correction, unchanged nuisance/contrast definitions and existing
UGR recovery/Trust helpers. Shell syntax checks passed.

A real FSL 6.0.7.19 `feat_model` probe on a synthetic UGR model produced two
180 × 13 designs at TR 1.615 and 1.710. Task columns changed (maximum absolute
change 1.2324606); contrast matrices and nuisance columns were numerically
identical. Probe output: `/private/tmp/rf1-TR-feat-model-2dm92rq2/result.json`.
The actual RF1 preprocessing and FEAT correction execute on linux2; no result is
claimed before that run and its QC complete.
