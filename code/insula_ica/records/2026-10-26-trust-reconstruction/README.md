# Trust reconstruction evidence, 2026-10-26

The completed sub-10402 Trust run 1 used exactly the historical FEAT nuisance
values. Its complete design and contrast matrices also match. Its BOLD midpoint
frame differs before FEAT masking, smoothing, intensity normalization, or FILM.
The reconstruction has not reproduced the historical statistical maps;
production remains stopped. Once validated, production should use eight
concurrent subjects.

## Evidence sources

- [Supplied linux2 diagnostic excerpt](linux2-preprocessing-excerpt.txt),
  preserved verbatim, including apparent transcription errors in a few path
  strings. Those strings are not independently validated filenames.
- [Structured summary](findings.json), manually transcribed from that excerpt
  and the earlier linux2 results supplied in this conversation. It identifies
  which evidence came from which source. Missing values remain unknown.
- Analysis code at commit `790cf515`: `inspect_trust_preprocessing.py`,
  `diagnose_film_replay.py`, and `rebuild_trust.py`, in `code/insula_ica/`.

These are user-executed linux2 results, not computations on the development Mac.
The full JSON remains on linux2 and was not supplied with the excerpt:

```text
/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/diagnostics/trust-preprocessing-evidence-20261027T032509Z-e54481.json
```

Its filename uses UTC (2026-10-27 03:25:09); the local date was October 26 EDT.
The corresponding terminal output was saved as
`/tmp/rf1-trust-preprocessing-evidence.bIMGoc`. The committed excerpt omits the
version section and sub-10369 results; neither is inferred from its omission.
The excerpt's SHA256 is recorded in `findings.json`.

## Exact nuisance identity

The historical FEAT directory retains `confoundevs.txt`, the actual nuisance
matrix used by FEAT. Recovery preserves a copy in `original_confounds/` and
points the new FSF to that copy. It does not substitute newly generated
fMRIPrep confounds or newly fitted TEDANA components.

| Run | Saved copy vs historical FEAT | Saved copy vs rebuilt FEAT | Frozen saved-file hash |
|---|---|---|---|
| sub-10402 Trust 1 | Exactly equal, 280 × 25, max difference 0 | Exactly equal, 280 × 25, max difference 0 | Matches |
| sub-10402 Trust 2 | Exactly equal, 280 × 26, max difference 0 | Rebuilt FEAT file unavailable | Matches |

The diagnostic uses `np.array_equal` after loading numeric values. For run 1,
this establishes element-for-element identity in the original and rebuilt
FEAT nuisance inputs, not merely the same column count or approximate agreement.
The saved-copy hash separately establishes unchanged bytes relative to the
recovery snapshot. The excerpt does not establish byte-for-byte equality of
the two FEAT text files, which may differ in whitespace, or independently
authenticate the missing external `confounds_tedana-24` file. Recovery uses a
copy, not the same filesystem object or original external pathname.

The earlier preflight reported preserved combined matrices for 241/241 eligible
Trust runs. That is not evidence that 241 reconstructed FEAT fits have used
them: only run 1 above has a completed rebuilt-matrix comparison in this excerpt.
Run 2 confirms preservation of its historical matrix but has no rebuilt design
yet. The sub-10369 exact nuisance comparison is not included in this excerpt.

**Standalone TEDANA is not an unresolved recovery dependency when the historical
combined nuisance matrix is reused.** Missing TEDANA metrics/mixing files do
not invalidate that reuse. The matrix preserves the original nuisance model,
including whatever TEDANA regressors that model contained. No new component
fit or extra regression is needed to recover those values. This does not prove
the original denoising choices were optimal; it establishes model reproduction.

## The discrepancy precedes FEAT prestats

For sub-10402 Trust run 1:

- No FSF setting differs after excluding the three relocated paths: BOLD,
  nuisance copy, and output directory.
- The full `design.mat` is 280 × 35 and `design.con` is 18 × 35; both comparisons
  report maximum absolute difference **0.0**.
- Both logs extract the same midpoint frame with
  `fslroi prefiltered_func_data example_func 140 1`.
- The saved settings specify 280 volumes, no volume deletion, no alternate
  reference, and no FEAT motion, slice-timing, or unwarping operation.

| Checkpoint | Spatial r | Relative RMS difference | RMS after best global scale | Common voxels |
|---|---:|---:|---:|---:|
| `example_func` | 0.995130 | 5.1591% | 5.1217% | 99,505 |
| `mean_func` | 0.997528 | 3.4584% | 3.4441% | 99,505 |

Both pairs share the same grid. Relative RMS uses the norm of the historical
image as the denominator, within the common masks. `example_func` is extracted
from the copied input before BET/smoothing/normalization in both logged command
sequences. Therefore at least part of the disagreement is already in the BOLD
entering FEAT. The global-scale fit does not remove it. This is an inference
from the retained images and command order; the deleted original 4D BOLD itself
cannot be compared directly. `mean_func` is overwritten after prestats, so it
is a later checkpoint.

Both logs contain 55 commands. Their SUSAN spatial sigma is the same
(`2.1231422505307855`, corresponding to the specified 5 mm smoothing). Derived
values differ:

| Quantity | Historical | Rebuilt |
|---|---:|---:|
| Intensity mask threshold | 656.3429688000001 | 657.5821289 |
| SUSAN brightness threshold | 2895.537048 | 2923.66735875 |
| Global intensity multiplier | 2.5901930714995984 | 2.5652713115785475 |

Those are data-derived quantities; their differences are consistent with the
different BOLD values and do not show a changed requested smoothing kernel.
The command diff also shows the expected relocated input pathname, because the
normalizer did not replace the extensionless BOLD names in FSL commands.

## Earlier linux2 results

These values were supplied in earlier terminal excerpts, rather than in the
new attachment. Raw copies of those earlier reports are not present locally.

| Trust run 1 | Maps passing comparison | Mask Dice | COPE10 r / RMS difference | ZSTAT10 r / RMS difference |
|---|---|---:|---|---|
| sub-10369 | 0/36; none missing | 0.994699 | 0.985898 / 16.633% | 0.965476 / 25.893% |
| sub-10402 | 0/36; none missing | 0.993877 | 0.987979 / 15.472% | 0.976189 / 21.662% |

The same-input UGR FILM control subsequently passed **34/34** model-defined
maps under **FSL 6.0.7.18**, with task-contrast relative RMS and maximum absolute
differences **0.0**, unchanged historical inputs, and no reported errors.
The terminal log was `/tmp/rf1-film-control.931ycY`. It supports reproducibility
of that retained-data statistical fit; it does not establish equivalence of
all Trust preprocessing.

The earlier input inventories reported:

- Fixed N=132 matches the corrected `include_final` cohort exactly.
- Trust: 241 accepted runs; UGR: 249 accepted runs.
- All 249 UGR runs retain a filtered-data FILM replay route.
- All 241 Trust runs lack source BOLD and filtered data at the checked original
  paths, but retain historical combined nuisance matrices.
- All 241 Trust runs have four raw magnitude echoes and JSONs, and usable
  Warpkit fieldmap/magnitude products and associations in BIDS; all 132
  participants have session-01 T1w files. No Warpkit rerun is currently required.
- Two pilot participants completed fMRIPrep. The Trust comparison failed before
  all three selected FEAT runs completed. No eight-subject production run started.

## Version certainty and next investigation

The actual rerun container reports **fMRIPrep 24.1.1** and its image hash is
frozen in the rebuild root. The archived cluster script
[`a19c46ef:code/fmriprep24-hpc.sh`](https://github.com/DVS-Lab/rf1-sra-linux2/blob/a19c46ef/code/fmriprep24-hpc.sh)
also specifies 24.1.1. That is evidence for the selected recovery configuration,
not proof of the original run's exact runtime. The attached excerpt contains no
original execution version record. Do not claim historical minor-version
identity solely from `fmriprep-24` directory names or the archived script.

The remaining investigation concerns BOLD generation upstream of FEAT: original
run-specific execution metadata, raw inputs and sidecars, fieldmap associations,
anatomical/template inputs, and multi-echo combination/runtime settings. Compare
the retained provenance at the already known RF1 paths before choosing another
rerun. Current evidence does not select one of those causes. No further TEDANA
component provenance request is needed to establish the retained nuisance fit.
Do not relax map agreement or launch the eight-subject expansion to bypass this
unresolved BOLD discrepancy. Existing completed work remains reusable.
