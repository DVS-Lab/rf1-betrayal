# Proceed with new 24.1.1 analysis inputs

Decision: 2026-10-27, following the user's acceptance of regenerated BOLD for
the new insula analysis despite nonidentical historical results. Stop the
container-version/reproduction investigation. Continue with fMRIPrep 24.1.1,
the fixed N=132 cohort, and all 241 originally accepted Trust runs.

This produces a new analysis derivative. It is not labeled an exact recovery
of the manuscript BOLD or statistics. The original model/confounds are retained.

## Execution

`start_trust_regeneration.sh` launches `regenerate_trust.py` quietly in the
background. Up to eight subjects run concurrently, each with 14 total processes,
eight OpenMP threads and a 28,000 MB scheduling cap. New runs use the recorded
master seed 20261027 and fixed skull-strip seed. This retains normal threaded
throughput; it is not the single-OpenMP-thread exact-repeat control recipe.

Completed fMRIPrep outputs for sub-10369 and sub-10402 are reused after grid,
time-count and TR checks. For other participants, a prespecified completed
24.1.1 control A output may be adopted after checking its image content hash,
runtime container hash, and reference geometry. The available additional subject
is sub-10478. Its subject directory is linked into the isolated rebuild root,
preserving the large control files in place. No output is selected by which
variant best matches historical statistics. Existing settings for reused runs
remain part of their provenance; those participants are not rerun solely to
standardize thread/seed settings.

The existing rebuild root and byte-preserved confound copies are retained.
Remaining work uses separate production scratch paths. A completed preprocessing
marker allows resuming without rerunning fMRIPrep. A partial FEAT directory is
preserved for inspection, never deleted/overwritten by an old wrapper. The
existing rebuild execution lock prevents concurrent runs of the old/new drivers.

## What remains required

- The actual combined nuisance values must match the historical FEAT matrix
  exactly, and the preserved copy's SHA256 must match its snapshot. No new
  fMRIPrep nuisance columns or newly fitted TEDANA components are substituted.
- The full task+nuisance design and contrast matrices must match. Only the
  BOLD path, preserved-confound path and FEAT output path change in the FSF.
- Spatial grid, TR, sample count, required output presence, and finite residual
  data are checked. FEAT fits the original joint model once; no sequential
  nuisance model or extra temporal filtering is added.
- Left, right and bilateral anatomical insula coverage is reported, using the
  same HarvardOxford atlas code as the ICA workflow. Constant residual voxels
  are omitted from functional coverage. Coverage below the existing 80% target
  is flagged for review.

Historical COPE/ZSTAT differences remain in their reports but do not stop
processing. The old exact-map pilot marker is neither required nor fabricated
by this production path. Historical-reference spatial correlation below 0.95
is a heuristic flag for visual registration review, not an exclusion criterion
or a replacement exact-reproduction gate. A flagged subject does not cancel
other subjects. Critical preprocessing/model errors are recorded per participant
while the remaining queue continues.

No participants/runs are silently removed or automatically repaired. Registration
and coverage review remain necessary before using the prepared data for ICA,
particularly given the large difference in one ordinary-setting control.

## Records and subsequent analysis

Each launch creates `code/insula_ica/records/linux2-trust-regeneration-*` with the
decision, container identity, input registry, commands, participant status, QC
metrics and error reports. It commits and pushes only that record at completion,
including failures. The background log contains short progress messages and
paths to detailed per-subject logs; it is not dumped into the terminal.

This prepares Trust BOLD and original-model FEAT residual inputs. UGR retained-data
FILM recovery and the shared TG+UG ICA follow as subsequent stages. The separate
native FILM residual-scale validation is not silently waived by acceptance of
new BOLD. No ICA or downstream statistical results are launched by this driver.

Validation: five new checks passed for exact nuisance identity, acceptance of
nonidentical historical statistics with unchanged designs, rejection of changed
designs, completed-control adoption without overwriting source data, production
resources/settings, and finite-data/coverage QC. All 15 existing Trust rebuild
checks passed. Static undefined-name/import checks and launcher shell syntax
checks passed. Actual full-cohort execution occurs on linux2.
