# UGR anatomy reuse and automatic unsmoothed FEAT recovery

The sm-0 audit found 249 missing original UGR BOLD paths under
`/ZPOOL/data/projects/rf1-sra-linux2/derivatives/fmriprep-24`.
The existing UGR sm-5 filtered data remain usable for the already running ICA
streams, but cannot supply pre-smoothing data. All 241 Trust sm-0 FEATs completed.
The user authorized regeneration of UGR with reused anatomy, a gated pilot,
automatic expansion, high concurrency and GitHub publication.

## Linux2 launch

```bash
bash code/insula_ica/start_ugr_anatomy_reuse.sh
```

The launcher uses the existing `tedana-26.0.3` Python environment for orchestration
and `/ZPOOL/data/tools/fmriprep-24.1.1.simg` for preprocessing. The image checksum
must match the recorded Trust image. It does not run a new TEDANA or Warpkit fit.
Required four-echo raw UGR and associated Warpkit BIDS fieldmaps are audited at
their exact known locations. Missing inputs are recorded, never substituted.

The pilot selects the first accepted single-run participant and first accepted
two-run participant before inspecting results. Each run uses an isolated one-run
BIDS view and scratch/output directory, retaining natural matching of run-specific
SBRefs with the **unchanged** upstream stock filter. The anatomy-only cache comes
from the same participant's completed 24.1.1 Trust reconstruction. The installed
sMRIPrep graph must contain no heavy anatomical processing, and the actual run log
must confirm precomputed registrations. A request for new anatomical normalization
terminates that processing child.

The three pilot runs must have correct grid, time count and TR, finite variable
BOLD, verified reference selection, at least 80% anatomical insula coverage on each
side and bilaterally, and historical reference spatial correlation >=0.90. Their
sm-0 FEATs must preserve the original complete design, contrasts and exact saved
TEDANA/nuisance matrix, have finite variable residuals and sufficient degrees of
freedom, and retain the residual/PE/autocorrelation intermediates. A pilot failure
stops expansion and publishes its record; these checks do not exclude subjects.
The 0.95 historical-frame threshold remains a visual-review flag. Historical
voxelwise equivalence of new preprocessing is not required.

After the pilot, functional preprocessing and FEAT overlap. Maximum queues are
32 subject workers and 64 FEAT workers. Each fMRIPrep invocation requests 8 total
threads, 4 OpenMP threads and 12,000 MB. A shared gate reserves 8 CPUs/12 GiB per
fMRIPrep job and 1 CPU/3 GiB per FEAT, with 12 GiB available-memory headroom.
CPU affinity and available memory/cgroup limits can lower simultaneous execution.
Completion time depends on the host and data; no overnight completion is implied
by the queue sizes.

UGR outputs go under `derivatives/insula_ica/ugr-anatomy-reuse-*`.
Regenerated BOLD is passed to sm-0 FEAT through an explicit identity/provenance
registry, never by replacing historical paths. A new generated UGR FSF is saved
beside the original missing-input FSF in the existing separate sm-0 root. Only
the input BOLD path, confound/output locations and requested smoothing setting
are changed; original task EVs/model and nuisance values are retained. Existing
Trust sm-0 completion markers/QC are validated and reused, not recomputed.

Records publish after the successful pilot, after the expanded recovery, and at
final combined FEAT completion or failure. Individual expanded-run failures are
preserved while other accepted runs continue. A combined normal sm-0 record
revalidates the full 490-run registry. Scientific QC review remains pending; this
launcher does not start ICA or change the current sm-5 ICA inputs. Publication
failures are recorded locally and do not stop ongoing processing.

Running the launcher again reuses the latest recovery root. Both the recovery and
sm-0 execution locks must be free, and lingering processing children block a new
submitter. Completed outputs are reused. Existing partial FEATs are preserved for
inspection; they are never deleted or automatically overwritten.

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/insula_ica/status_ugr_recovery.py --brief
```

## Development validation

63 selected local tests passed, including the full mocked 132-participant /
249-run automatic expansion, failure before expansion, unchanged UGR stock-filter
views, UGR-specific Warpkit association checks, real NIfTI header/finiteness QC,
explicit regenerated-BOLD provenance, preservation of original FSFs and completed
Trust outputs, and regression coverage for existing Trust/ICA helpers. The shell
launcher passes `bash -n`. These checks do not imply a completed RF1 UGR pilot;
the actual 24.1.1 pilot and numerical FEAT checks execute on linux2.
