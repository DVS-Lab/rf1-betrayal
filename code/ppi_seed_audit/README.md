# Original PPI seed provenance audit

This is a read-only audit of the original project, not a PPI or FEAT rerun.
It distinguishes a seed label in a filename from evidence about the input
actually used by a saved analysis. Existing scripts, masks, physiological files,
models and group outputs are never modified.

## Linux2 launch

```bash
bash code/ppi_seed_audit/start.sh
```

The launcher uses the existing `tedana-26.0.3` Python, single-threaded numerical
libraries and `nohup`. Results publish automatically to `insula-ica` in
`code/ppi_seed_audit/records/linux2-PPI-seed-audit-<timestamp>/`. It stages only its
own record. A lock prevents duplicate running audits. Rerunning after completion
creates another inventory; it does not restart or change production analyses.

```bash
/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python \
  code/ppi_seed_audit/audit.py --brief
```

## Evidence collected

- Relevant current scripts/templates and up to 100 Git commits affecting those
  scripts, with original content snapshots and line-numbered seed/extraction
  excerpts. Scripts are inspected, never executed.
- Saved N132 group designs under the original `derivatives/fsl/L3*` directories,
  following insula-like names or referenced inputs through L3 → L2 → L1. Group
  folder names do not have to say PPI: naming mistakes must not hide the chain.
  Unrelated N132 groups are indexed but not followed. Original insula-like
  participant directories are inventoried even if no surviving group references
  them, with group-chain membership marked separately.
- Original L1 `design.fsf`, matrices/contrasts and physiological custom-EV paths.
  The cohort registry is the already verified N132/490-run manuscript registry.
  Existing excluded-run directories within those participants are marked as
  such; no excluded run is added to an analysis or numerical BOLD check.
- Raw physiological snapshots in `custom_timing_files/evN.txt`, the external
  physiological file referenced by the FSF, and disagreement between them. FSL
  normally copies custom timing inputs there for reference. The snapshot takes
  priority when present, because an external file may have been overwritten.
- Agreement with the physiological column actually stored in `design.mat`.
  For supported direct, unconvolved and unfiltered EVs, this requires agreement
  up to positive centering/scaling with relative residual below 1e-6. If direct
  text comparison fails, an explicit float32 input representation is checked:
  native `feat_model` reads shape-2 custom values using `fscanf("%f", &tmpf)`
  with `float tmpf`. Both comparisons are retained; the cutoff is unchanged.
  This does not relax the separate raw mask-mean identity check. Unsupported
  filtering/basis expansion/orthogonalization is explicitly left unresolved.
  A matching series alone does not identify a spatial mask.
- Bounded original FEAT and project wrapper-log excerpts, parsed `fslmeants`
  commands, current mask candidates and recovered historical mask versions.
  Literal commands are distinguished from unresolved shell-variable commands.
  Source code text alone is not treated as a record of execution.

Scopes are the original project directories and exact saved references. Stale
same-project FSF/physiological paths may be relocated by their project-relative
suffix for read-only inspection; the original path and relocation are recorded.
There is no filesystem-wide search and no substitution of new-version or rebuilt
BOLD for missing original inputs.

## Numerical candidate checks, when original BOLD survives

By default, at most four distinct **exact FSF-referenced original BOLD images**
are read serially. Sampling covers task/seed-label families before any numerical
fit is inspected. One 4D read supplies all candidate-mask means for every saved
series referencing that image. Available RAM is checked first on Linux.

Only binary/integer-label masks sharing the exact voxel lattice are compared.
An exact field-of-view embedding may use the same physical voxel centers on a
larger image; no interpolation or resampling is performed. Means use float64
accumulation from float32 image data. The implementation was checked against
native FSL 6.0.7.19 `fslmeants` on a synthetic integer-label mask; maximum
difference was 4.85e-6 (see `native-mean-validation.json`).

Raw agreement requires `allclose(rtol=1e-7, atol=1e-5)` and centered RMS error
below 1e-4 of the observed series SD. Correlation or affine agreement alone is
not called a raw match. Candidate support for the saved model additionally
requires agreement between the selected raw series and the saved design column.
Equivalent mask supports retain all filename aliases; numerical agreement cannot
distinguish different filenames holding the same mask. Missing original BOLD or
absent contemporaneous extraction evidence can leave the exact seed unresolved.

`--max-bolds 0` disables these optional image reads; the script/design/timecourse
audit still runs. No BOLD image is uploaded. Candidate mean timecourses, input
identities and comparisons are saved for review.

## Outputs

- `state.json`, `report.md`: phase/counts, missing sources and summary.
- `analysis_chain.json`, `saved_designs.json`: actual saved links and group usage.
- `physiological_EVs.json`, `physiological_index.tsv`: external/cache paths,
  model membership, accepted-run flags and whether each case feeds a group.
- `script_history.json`, `script_evidence.json`, `log_evidence.json`,
  `extraction_commands.json`: provenance and extraction-command evidence.
- `candidate_masks.json`, `candidate_timecourse_comparisons.json` and
  `numerical_checks/`: spatial candidates, raw/affine matches and reconstructed
  one-dimensional series where possible.
- `identical_series_across_seed_labels.json`: different labels sharing exactly
  the same numerical physiological series within a run.
- `source_manifest.json` plus split `source-files-*.tar.gz`: exact source bytes
  identified by SHA256. Archive members are `files/<sha256>`; the manifest maps
  them to original filenames and Git commits. Splitting at 32 MiB uncompressed
  keeps individual artifacts below GitHub file limits and avoids case-only
  filename collisions. Files above the 8 MiB small-source limit are reported,
  not silently copied. Only small text/model/timecourse sources and 3D masks are
  archived, never participant 4D data.

Local validation: six targeted tests pass, including an end-to-end original
L3/L2/L1 fixture, actual project templates, detection of an overwritten external
series, recovery from FEAT's preserved EV snapshot, numerical mask attribution,
exact-grid embedding and content-addressed archive integrity. Native FSL mean
agreement, shell syntax and CLI help were also checked. These are fixture tests,
not claims about Jacob's exact seed mask. See the dated audit record's `review/`
for the actual saved-design comparison and remaining provenance limits.
