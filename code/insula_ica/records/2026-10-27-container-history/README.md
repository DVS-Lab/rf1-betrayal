# fMRIPrep container history and repeatability controls

Reviewed 2026-10-27, America/New_York. The table dates committed script
references, not the first actual execution or creation of a container file.
Exact script/config snapshots and their SHA256 values are saved in
[references.json](references.json) and `sources/`.

## What DVS-Lab Git history establishes

| Image filename | Dated evidence | Meaning |
|---|---|---|
| `fmriprep-24.1.1.simg` | [multiecho-pilot, 2024-12-18](https://github.com/DVS-Lab/multiecho-pilot/commit/5033f79d7085cfc30b3f7012d72f65dc466b3809) | Linux2 wrapper selected this image. |
| `fmriprep-24.1.1.simg` | [rf1-sra-data, 2025-03-14](https://github.com/DVS-Lab/rf1-sra-data/commit/cc81807f790284661c9added6615178ed039946f) | RF1 cluster wrapper selected 24.1.1, 28 total threads, Warpkit fieldmaps prepared separately, no fixed-seed flags. |
| `fmriprep-24.1.1.simg` | [rf1-sra-linux2, 2025-11-13](https://github.com/DVS-Lab/rf1-sra-linux2/commit/a19c46efdd89ea80f50fc12f3e891455e81cffc8) | Session-filtered cluster wrapper selected 24.1.1, 14 total threads, output `derivatives/fmriprep-24`, no fixed-seed flags. |
| `fmriprep-25.2.4.simg` | [rf1-sra-linux2, 2026-04-14](https://github.com/DVS-Lab/rf1-sra-linux2/commit/1ccfc94c75aaa0844a72474cc130b49b4dbeffaa) | Linux2 wrapper explicitly selected 25.2.4. |
| `fmriprep-25.2.5.simg` | [rf1-sra-linux2, 2026-04-14](https://github.com/DVS-Lab/rf1-sra-linux2/commit/35e6cf0dc8cb41f3d299bc2349c8e39e887bea19) | Explicit image replacement from 25.2.4 to 25.2.5. |
| `fmriprep-warpkit.sif` | [rf1-sra-linux2, 2026-05-04](https://github.com/DVS-Lab/rf1-sra-linux2/commit/3714c3a3508cedbce8997796959c0584c6aa05f2) | Test wrapper used `--me-use-warpkit` and wrote `fmriprep-warpkit-testing`; comment identifies the DVSneuro `warpkit-medic` fork. |
| `fmriprep-warpkit-andrew.sif` | [rf1-sra-linux2, 2026-05-16](https://github.com/DVS-Lab/rf1-sra-linux2/commit/e17334e7114db0047d660a40f20a3464ca5cb97a) | Wrapper changed the test image to Andrew's version, still with integrated Warpkit and separate validation outputs. |

The [official changelog](https://fmriprep.org/en/25.2.5/changes.html) dates
25.2.4 to January 14, 2026 and 25.2.5 to March 10, 2026. The RF1 wrapper's
explicit upgrade is in April. These findings strengthen 24.1.1 as the
historically supported candidate, while the newer images can serve as
comparison arms. Closer agreement from a newer image would not prove it was
used for the manuscript.

The integrated Warpkit images are development/test pipelines. Their filenames
do not establish their underlying fMRIPrep version or exact fork commit.
`run_bold_controls.py` inventories all five existing linux2 images, recording
content SHA256, runtime package versions and container labels. It runs only
the requested 24.1.1/25.2.5 comparisons; it does not refit Warpkit or use the
development images as historical substitutes.

The archived scripts are configuration evidence, not subject/run execution
proof. Historical image hashes and original run-specific configs remain
unavailable. Current file mtime is not a reliable build/first-use date.

Search coverage: indexed GitHub filename searches for the five image names;
all locally available Git refs for `rf1-sra-linux2`, `rf1-sra`, `rf1-sra-trust`,
`rf1-sra-ugr`, `rf1-sra-sharedreward`, `multiecho-pilot`, and `rf1-betrayal`;
and the older `rf1-sra-data` cluster-wrapper path history and exact dated file
via the GitHub connector. See [searched-histories.json](searched-histories.json).
No matching image-name changes in some repositories means only that those
searched paths/history did not contain them. It does not establish non-use.

## Should independent fMRIPrep runs be identical?

Exact voxel equality is not an established expectation under the earlier
uncontrolled settings. Both the
[24.1.1 documentation](https://fmriprep.org/en/24.1.1/usage.html#specific-options-for-ants-registrations)
and [25.2.5 documentation](https://fmriprep.org/en/25.2.5/usage.html#specific-options-for-ants-registrations)
describe repeatability controls: `--random-seed`, `--skull-strip-fixed-seed`,
and `--omp-nthreads 1`. Our previous pilot saved varying random seeds,
`skull_strip_fixed_seed=false`, and eight OpenMP threads. The archived commands
did not request fixed seeds either. The strict historical map-equality check
is useful as a diagnostic, but was too strong as the sole acceptance rule for
that configuration. A failed exact comparison alone does not establish a
wrong preprocessing version. The observed 4–8% frame differences still need
an empirical repeat-run baseline before being attributed or dismissed.

Even with repeatability controls, compare numerical voxel values rather than
compressed file bytes or report/header timestamps. Matching software version
labels alone do not control image contents, library dependencies, templates,
inputs, seeds, or thread settings. Cross-version equality is not presumed.

## Implemented experiment

Select the existing sub-10402 pilot plus the first additional QC-eligible,
Warpkit-ready subject with both accepted Trust runs (normally sub-10478).
Selection is fixed before examining new outputs. Eight independent subject
jobs use fresh, distinct output/work directories:

1. Sub-10402: two ordinary 24.1.1 runs with eight OpenMP threads and no fixed seed,
   compared with each other and its retained pilot BOLD from the same image SHA256.
2. Both subjects: two independent 24.1.1 controls with the same master seed,
   fixed skull-strip seed, and one OpenMP/ITK/BLAS thread.
3. Both subjects: one 25.2.5 control with those same requested settings.

All arms retain Trust-only selection to isolate the repeat/version questions.
Restoring all-task selection would be a separate experiment, not a concurrent
change to the version arm. The retained BIDS Warpkit fieldmaps are reused.
The BIDS input and shared template cache are mounted read-only. Raw/anatomical
images, relevant metadata/fieldmaps, existing pilot BOLD and reference frames,
templates and tested image contents are hashed before/after. The image-content
identity check prevents an unnoticed replacement of the old pilot container.
Saved execution configs verify the requested seed/thread controls were applied.

Comparisons include every time point for the repeat/version pairs: raw relative
RMS, RMS after global scaling, voxel-demeaned relative RMS, and voxelwise temporal
correlation percentiles. Each arm's midpoint frame is compared with the retained
historical FEAT reference. Different grids are reported without resampling.
The exact same subject/run paths under the versionless `derivatives/fmriprep`
root are also inspected; the April 2026 wrapper explicitly used that root.

This is not an exact reconstruction of unknown historical seeds or an exhaustive
estimate of stochastic variability. No FEAT fit, standalone TEDANA, new Warpkit
fit, passing pilot marker, or production expansion is launched. The experiment
reports differences and does not automatically select a preferred version.

The quiet background launcher writes small progress messages to a `/tmp` log,
saves full results/configs/case logs under `code/insula_ica/records/`, and commits
and pushes only that record on completion, including failures. Large images/work
files stay under ignored derivatives/scratch roots. Git publication failure is
reported in the log and retained record rather than losing the analysis results.

Validation: seven tests passed, including full-series precision/temporal checks,
fresh output/work paths, and an end-to-end eight-case fixture with mocked
preprocessing. Original anchor hashes remained unchanged. No real fMRIPrep
container was executed on the development Mac; the actual controls run on linux2.
