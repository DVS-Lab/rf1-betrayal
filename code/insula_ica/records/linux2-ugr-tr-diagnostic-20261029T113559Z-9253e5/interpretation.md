# Timing discrepancy in sub-10817 UGR run 1

The published `evidence.json` establishes:

| Source | TR (seconds) | Volumes |
| --- | ---: | ---: |
| All four current raw BIDS echo JSONs | 1.710 | — |
| All four current raw NIfTI headers | 1.710000038 | 240 |
| Original FEAT `design.fsf` | 1.615000 | — |
| Retained sm-5 `filtered_func_data` header | 1.615000010 | 240 |

The saved nuisance matrix hash still matches the accepted registry. The historical
pre-FEAT BOLD and its JSON are absent at their recorded paths. This is an actual
0.095-second discrepancy, approximately 5.88% relative to the modeled TR, not
floating-point rounding or a milliseconds/seconds conversion.

The retained FEAT header does not independently establish the acquisition TR.
Installed local FSL 6.0.7.19 `tcl/featlib.tcl` lines 5480–5486 reads `pixdim4`
from `filtered_func_data`; if it differs from `fmri(tr)` by more than 0.01, it
rewrites the header `dt` using the FSF value via `fslhd`/`fslcreatehd`.
Linux2's recorded FSL is 6.0.7.18; its exact installed source was not fetched in
this diagnostic. The verified local behavior is sufficient reason not to treat
the post-FEAT header as independent proof of original scan timing, but does not
prove that this historical run underwent that operation.

Current raw metadata consistently support 1.710 seconds. They do not alone prove
that these are unchanged historical inputs, nor establish whether the original
model's 1.615 seconds was wrong. Original DICOM acquisition metadata, retained
conversion/preprocessing provenance, or task scanner-trigger timing should resolve
that question. Correct task event alignment must be checked if the modeled TR is
changed. Do not change raw headers/JSONs to 1.615 merely to pass the audit, and do
not claim that a model at 1.710 preserves the historical temporal design.

Status remains 489/490 completed unsmoothed FEATs. This one run never entered the
new UGR fMRIPrep workflow. Its already accepted historical sm-5 residuals were used
in the earlier ICA streams, so any established timing/model correction would also
need to be tracked for those inputs. No source metadata, models, exclusions or
analysis outputs were changed by this interpretation.
