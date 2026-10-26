# Exploratory insula ICA for the betrayal revision

This self-contained analysis uses **left, right and bilateral anatomical insula**,
each with one shared TG+UG decomposition. Five components per ROI are primary;
three and seven are spatial-stability checks. Component numbers do not correspond
between ROIs or dimensionalities. Existing FEAT models, results, upstream
preprocessing and the N=132 participant list are never edited.

## Current status and input gate

The implementation can be tested locally, but **production is blocked until the
linux2 audit passes**. No participant ICA/connectivity/task results have been
computed during implementation. `/ZPOOL` is not available on the development Mac.

Source inspection on 2026-10-26:

- Betrayal main: `affb3b30ad9b911232a99e32df6459403c84d4ea` (local clone originally
  `6aba28a25c7256b2bcca03339f3996da9fa200f3`). Read root and code READMEs.
- Upstream local clone: `72a5a160a8636cff770ea93116cb13bab64cf609`; remote main:
  `02faf289fcf04908c4a57bdc49ddf9e3b7498c97`. Current confound generator inspected
  at that remote commit; it is **not** evidence that legacy `-24` files match it.
- Archive: `8d127f18d5930803410e716a4b62aac7d4aa4beb`; inspected
  [`leech12_dual_regression`](https://github.com/DVS-Lab/duke-archive/blob/8d127f18d5930803410e716a4b62aac7d4aa4beb/leech12_dual_regression)
  and [`leech12_dual_regression_conf`](https://github.com/DVS-Lab/duke-archive/blob/8d127f18d5930803410e716a4b62aac7d4aa4beb/leech12_dual_regression_conf).
  The latter appends confounds to stage 2. This implementation uses explicit
  matrix operations and keyed run records rather than its parallel positional
  lists and `paste` command.
- Both task templates specify TR=1.615 s, `ndelete=0`, and disabled FEAT temporal
  high-pass filtering. The actual NIfTI TR and saved FSF must agree at audit.
- Local FSL 6.0.7.19, MELODIC 3.15; local help inspected for MELODIC and fsl_glm.
  Linux2 versions/help are captured again by the audit.
- `code/sublist_n132.txt` contains 132 unique participant IDs. The saved
  `master_subject_table_full.csv` was absent from the local checkouts, current
  GitHub main tree and available betrayal Git history. It may be on linux2.
  **Do not regenerate its outlier thresholds on a different cohort or infer run
  inclusion from image availability.** Configure its existing saved location.
- Independent Harvard–Oxford cortical maxprob thr25 mask, XML Insular Cortex
  index 1 / NIfTI label 2. Native 2 mm counts: left 1,160, right 1,181,
  bilateral 2,341. These are **not** final BOLD-grid coverage counts.

The script-derived source candidates, not yet live-data verified, are:

```text
/ZPOOL/data/projects/rf1-sra-linux2/derivatives/fmriprep-24/sub-<id>/ses-01/func/
  sub-<id>_ses-01_task-<trust|ugr>_run-<1|2>_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
/ZPOOL/data/projects/rf1-sra-linux2/derivatives/fsl/confounds_tedana-24/sub-<id>/
  sub-<id>_ses-01_task-<trust|ugr>_run-<1|2>_desc-TedanaPlusConfounds.tsv
```

Their matching fMRIPrep confounds, BOLD JSON and MNI brain masks must exist.
The audit checks those **exact** candidates against each original activation
FEAT `design.fsf`, including its BOLD and confound paths. It never substitutes
copes, zstats, FEAT residuals, new upstream derivatives or TEDANA denoised files.
An optimally combined fMRIPrep image and a TEDANA-denoised image are distinct:
`part-mag` / `desc-preproc` alone does not establish either provenance.

## Denoising decision

The primary plan is **joint nuisance cleaning before ICA**, followed by
confound- and task-adjusted dual regression. Adding confounds only after ICA
cannot recover a useful component basis if noise dominated discovery.
[Lindquist et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6865661/)
warn about sequential nonorthogonal projections, not nuisance regression before
ICA in itself. We use one nuisance subspace Q and consistently project both
BOLD and subsequent task columns with R = I − QQᵀ.

The exact existing legacy confound matrix is used, with an intercept. Its prefix
must match the inspected upstream generator layout: six aCompCor columns, six
motion terms, FD, existing cosine terms and NSS spikes. Its remaining columns
are preserved, including TEDANA noise terms. Drift is included **in this same
projection**: existing cosines are retained; only when none exist is a 128 s DCT
basis added. No sequential second high-pass, event regression for discovery,
extra smoothing, new global-signal regression, new censoring threshold, or TEDANA
rerun is performed. NSS spikes remain on the original time grid; initial NSS
rows are omitted only from ICA discovery. Component task/connectivity models
retain the original event timing and volume count.

The config's provenance review must establish what generated the selected BOLD,
which echo/optimal-combination operation was used, whether it was already
TE-denoised, and the source of the appended TEDANA noise columns. If those noise
terms were already removed, **stop**: this primary implementation is specifically
for the legacy preprocessed-BOLD + confound-regression route. Do not set its
review flag merely to make a gate pass. It does not authorize replacing data or
removing columns from already denoised data without a revised plan.

## Run on linux2

Python >=3.10 with NumPy, SciPy, NiBabel and Matplotlib is needed, together with
FSL `melodic`, `feat_model`, `fsl_glm` and the standard atlases. Install
`requirements.txt` in a suitable existing environment if needed. Set `PYTHON`
to that interpreter; shell `python3` can resolve to an unrelated installation.
Work is sequential and currently loads one whole 4D run; allow ample RAM
(roughly several times the uncompressed run size) and disk for cleaned images
and per-run connectivity maps. No scheduler assumptions or upstream changes.

From your linux2 betrayal checkout:

```bash
git fetch origin
git switch --track origin/insula-ica
cp code/insula_ica/config.example.json code/insula_ica/config.local.json
export PYTHON=/absolute/path/to/your/python
"$PYTHON" code/insula_ica/insula_ica.py --config code/insula_ica/config.local.json audit
```

The **first audit is useful with the review fields blank**. It reports blockers,
exact paths, all 528 potential cohort/task/run slots, original QC exclusions,
accepted counts, geometry/TR, nuisance rank, NSS count, residual dof, EV timing
and mask coverage. `audit/report.md` and `audit/manifest.tsv` are the first files
to inspect. `master_qc` can be an absolute path to the saved linux2 CSV; no need
to copy it into Git. Missing/unknown QC produces a blocker, never implicit
inclusion. All N=132 identities remain required; missing data stops production
and is disclosed for scientific review rather than silently reducing N.

Fill `provenance_review` from the original sidecars/commands and saved EV coding,
set a **fresh output_root**, then repeat `audit`. Fresh output roots preserve
failed attempts and prevent stale outputs from being mixed. Relative config
paths resolve against the betrayal repository, not your current directory.
The live upstream Git commit is recorded. Keep the legacy `-24` derivative
settings unless a separately reviewed provenance decision changes the source.

Configure independent MNI network templates before the final `summarize` stage:
provide one 3D map, name and DOI/version per entry in `network_templates`.
No atlas is downloaded or selected based on the observed task effects.
Templates must already be in the same MNI reference; only grid interpolation is
performed, not registration. Missing templates leave this stage blocked.

After a passing audit, run:

```bash
"$PYTHON" code/insula_ica/insula_ica.py --config code/insula_ica/config.local.json prepare
"$PYTHON" code/insula_ica/insula_ica.py --config code/insula_ica/config.local.json --scope pilot ica
"$PYTHON" code/insula_ica/insula_ica.py --config code/insula_ica/config.local.json --scope pilot regress
"$PYTHON" code/insula_ica/insula_ica.py --config code/insula_ica/config.local.json --scope pilot summarize
```

The deterministic pilot spans up to five participants across nuisance rank and
NSS burden and retains their original valid runs from both tasks. Review
`pilot/regression_qc.tsv`, task EV timing/amplitudes, reports for all nine ICA
solutions, mask loss, component timecourses, sign convention and template
correspondence. A successful synthetic pilot is **not** a real-data pilot.
The full ICA command requires completed real-data pilot regressions.

```bash
"$PYTHON" code/insula_ica/insula_ica.py --config code/insula_ica/config.local.json --scope full ica
"$PYTHON" code/insula_ica/insula_ica.py --config code/insula_ica/config.local.json --scope full regress
"$PYTHON" code/insula_ica/insula_ica.py --config code/insula_ica/config.local.json --scope full summarize
```

Input identity, lightweight file checksums, BOLD size/mtime, software versions,
config and pipeline source checksum are rechecked between stages. Full BOLD
checksums are deliberately not calculated for hundreds of private 4D images.
Commands, stdout/stderr, exit codes and source EV/FSF identities are preserved.
A changed analysis config or input requires a fresh audit/output root.
Independent network-template entries/notes are the only config exception: they
may be supplied later because they do not change ICA or task/connectivity fits;
the summary saves their exact identities and citations. Do not rerun a
partially written processing stage in place; preserve it and use a new output
root. Existing complete stages are not overwritten.

## Estimation details

1. Masks: nearest-neighbor resample the independent 2 mm HO anatomical label,
   split by MNI x<0 / x>0 (not array direction), intersect with the fixed common
   MNI brain-mask coverage across all accepted runs. This is anatomical cortical
   gray matter coverage, not an activation/statistical mask. Record final and
   lost voxels, retained fraction and lost-voxel MNI bounds. The coverage gate is
   80%; inspect the lost maps for anterior/ventral loss even above that threshold.
   Nonfinite/zero-variance data within whole-brain coverage stop preparation.
2. ICA: pooled TG+UG temporal concatenation separately for each ROI. Preserve
   every accepted steady-state frame and run. Build **one input per participant
   per task** by concatenating that participant/task's valid runs. Each run is
   mean-centered/variance-normalized voxelwise, then weighted by
   `1/sqrt(number_of_valid_subject_task_runs * steady_frames_in_this_run)`
   before concatenation. This equalizes run covariance contributions inside each
   participant/task input. Standard non-MIGP group MELODIC then performs its
   individual PCA whitening on each input: each participant/task supplies the
   same number of whitened PCs. This is session balancing, not a claim that
   unwhitened group covariance is exactly the weighted raw covariance.
   Disable MELODIC's voxel variance normalization with `--vn`, which the local
   help defines as **switch off**, so it does not undo within-session run weights.
   Disable MIGP to retain the standard individual-PCA/group-concatenation path.
   Disable automatic mask updating with `--update_mask` (the documented switch
   off), and require the resulting mask to equal the audited mask exactly.
   Save exact input/run boundaries, weights and all native MELODIC outputs/reports.
   ICA is stochastic; freeze saved solutions rather than expect bitwise agreement
   from an unsupported random-seed flag. A requested dimensionality that returns
   fewer maps is a blocker. Inspect convergence/restarts in the saved logs.
   Group-report spectra include concatenation boundaries; use saved **original-run**
   timecourses for temporal artifact review rather than treating these boundaries
   as continuous physiological time.
3. Sign: the voxel with largest absolute loading within each anatomical mask is
   positive (ties follow array order). Save original and oriented maps and the
   sign vector; no task coefficient or connectivity template determines sign.
   All 3/5/7 spatial pair correlations are reported; splits/merges need qualitative
   interpretation rather than forced one-to-one component matching.
4. Stage 1: simultaneously spatially regress all five maps plus a spatial
   intercept into the cleaned, percent-signal BOLD. Save the resulting five TCs
   and their scaling. TCs are standardized to one run-specific temporal SD.
5. Stage 2: regress all five TCs into whole-brain BOLD, controlling the same
   nuisance subspace and all existing task EVs. The implementation uses the
   Frisch–Waugh–Lovell equivalent of the full simultaneous model and saves the
   full design `[five TCs, projected task EVs, Q]`. Rank, condition and residual
   dof are checked. Maps are partial regression coefficients in percent signal
   per TC SD, **not bivariate correlations or voxelwise inferential z maps**.
   Collinear nuisance/empty task columns are represented by their column space;
   all five component effects and task contrasts must remain estimable.
6. Task effects: use original FSF/EV definitions rendered by FSL `feat_model` in
   analysis output only; the saved task-only FSF disables appended confounds,
   which are handled jointly by Q. No new HRF/event model is invented. TG is
   negative contrast 10 `rec-def`, retaining all partner conditions. UG is
   contrast 11 across all four pmods (partner × endowment), reversed only if
   the existing amplitudes increase with fairness. The sign review is required:
   current event-generation history includes renaming conventions, so filenames
   alone are insufficient. Store the original contrast magnitude convention
   (sum, not redefined average) and signed coefficients.
7. Run task uncertainty: OLS coefficients with Newey–West temporal covariance,
   lag `floor(4*(T/100)^(2/9))` (at least 1), with a residual-dof correction.
   This avoids treating fMRI samples as independent for reported within-run SEs;
   it is not a replacement for FEAT FILM or a claim of identical first-level
   uncertainty. Equal-run averages form participant effects; reported participant
   HAC SEs assume independent runs. Group t tests/95% CIs use empirical variation
   across participants, giving each participant one observation.
8. Group connectivity: equal-run means within participant/task, then equal-task
   participant means, then equal-participant group mean. Save participant/task,
   participant-pooled and group maps for every IC. Network correspondence is a
   signed spatial correlation outside the defining ROI with independent templates.
   Maps and correspondence are descriptive; there is no voxelwise inference or
   automatic assignment of a psychological label. Different maps alone do not
   constitute a formal component-difference test.
9. Multiplicity: **30 two-sided primary task tests** (3 ROIs × 5 ICs × 2 tasks),
   Bonferroni corrected as one family. Sensitivity dimensions receive no task
   discovery tests. No formal TG–UG interaction, cross-ROI component matching,
   component-wise gPPI or behavioral search is implemented.

## Outputs and interpretation

All generated outputs are ignored under `code/insula_ica/output/` by default.
If using another output location, keep private BOLD and bulky intermediates out
of Git. Select lightweight tables/figures deliberately for any later commit.

- `audit/`: manifest/exclusions, source/command provenance, task designs and EV
  timing/amplitudes, nuisance bases, geometry/coverage masks and dry-run report.
- `prepared/`: joint-cleaned 4D images and orthogonality diagnostics.
- `pilot/` / `full/`: 3/5/7 MELODIC solutions for each ROI, sign files, spatial
  stability matrices, discovery weights, all five component TCs and partial maps
  per run, signed run/participant effects and model diagnostics.
- Final summary: all 30 estimates/N/SE/CI/p/corrected-p, descriptive group maps,
  independent-template correspondence, component anatomy, two composite figures
  and `results_summary.md`. Component spectra/overlays are also in native
  MELODIC reports. Retain potential artifact components and inspect them.

This exploratory analysis can characterize heterogeneity in independent
anatomical insula and its task responses. The component basis is learned from
these task datasets: task effects are **not independent validation**, a unique
encoding claim, a replication of the previous aIns–vmPFC gPPI, or a test of its
survival with an independent seed. Null/ambiguous components must be shown.
TG and UG coefficients have different units and must not be directly subtracted
as a task interaction. A targeted vmPFC follow-up requires separate scientific
review; it is not implemented here.

Methodological sources:
[Leech 2012](https://doi.org/10.1523/JNEUROSCI.3689-11.2012),
[Leech 2011](https://doi.org/10.1523/JNEUROSCI.5626-10.2011),
[Lindquist 2019](https://doi.org/10.1002/hbm.24528),
[FSL dual regression](https://fsl.fmrib.ox.ac.uk/fsl/docs/resting_state/dualregression.html).
The supplied brief and complete archived scripts were available; full-text
retrieval of the two Leech papers was blocked by publisher/PMC access responses
in this session. Verify their full methods before manuscript methods are finalized.

## Validation

```bash
"$PYTHON" -m unittest discover -s code/insula_ica -p 'test_*.py' -v
"$PYTHON" code/insula_ica/test_insula_ica.py --integration /absolute/new/scratch/folder
```

Numerical tests check joint-projection task estimates, equivalence of partial
and simultaneous regression, contrast nonestimability, coordinate-based
hemisphere masks, sign orientation and the 30-test family. The integration test
creates synthetic positive 4D data, confounds, NSS, original-template EV designs,
three participants, both tasks, excluded second runs and fake reference maps.
It runs the actual installed `feat_model`, all nine real MELODIC solutions,
both regression stages and summaries for pilot and full scopes. It does not use
private imaging or establish scientific validity on RF1 data.

Development validation (2026-10-26): six numerical tests passed. The real-FSL
synthetic integration completed pilot and full scopes, all 18 MELODIC calls,
all regressions, 30-test summaries, template correspondence and both figures.
Stage-1 TCs also matched an independent `fsl_glm --demean` calculation within
1e-4 after the documented SD scaling. Fixed-mask equality was checked for every
solution. Figure layouts were visually inspected. These are software tests,
not RF1 scientific results.

For an audit-only Linux2 startup, `audit_linux2.py` locates an existing saved QC
CSV, verifies it against the unchanged N=132 list, creates a fresh timestamped
output/config, checks the active Python/FSL dependencies and runs only `audit`.
Use the upstream-documented interpreter
`/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python`.
An explicit file can be supplied as `--qc /existing/path/master_subject_table_full.csv`.
If supplied with `--qc-snapshot-stdin`, it reads a gzip/base64 attachment-derived
cohort-only CSV and searches for matching saved run flags; if none is found, it
restores that minimal snapshot in the ignored output directory. Neither the full
attachment nor the snapshot is published to Git. The full table can include more
`include_final` participants than N=132; the unchanged manuscript list always
determines analysis identity. This startup helper does not waive denoising or
UG-sign review gates and does not install packages or alter shared environments.
