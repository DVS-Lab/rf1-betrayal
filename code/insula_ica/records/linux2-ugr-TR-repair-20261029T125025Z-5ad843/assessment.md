# Assessment of completed timing correction

All three repair/control fits completed, and the combined unsmoothed registry
contains **490/490 completed FEATs**. The corrected case has finite residuals,
195 residual degrees of freedom, unchanged nuisance values and contrast
definitions, and a regenerated task design. Its left/right/bilateral insula
coverage is 100%, with no automated BOLD review flags. Across all 490 unsmoothed
FEAT QC reports the minimum per-run coverage is 100% left, 99.53% right and
99.76% bilateral. These are per-run coverage minima, not a common-mask intersection
calculation or completed visual registration review.

The audit of all 490 accepted runs found **one** raw-acquisition/model TR mismatch
and zero unresolved checks. Both sub-10817 Trust runs consistently have raw/model
TR 1.615 seconds. Only the accepted UGR run has current raw TR 1.710 versus its
historical model's 1.615. Its protocol/series name still contains `TR1615`, so a
protocol name is not sufficient evidence of the actual timing setting.

The same-BOLD sm-5 comparison establishes a large first-level effect of changing
the modeled TR. UGR contrast 11, offer (un)fairness, has COPE spatial r=-0.1253 and
ZSTAT r=-0.1180 between the wrong-TR and corrected models, with relative RMS
differences of 140.14% and 143.93%. This is not a preprocessing-version comparison:
the BOLD, smoothing, event files and nuisance values were held constant. Over a
240-volume acquisition, the assumed duration differs by approximately 23 seconds.
Weak negative spatial correlation does not establish a uniform sign reversal.
These single-run map metrics do not establish the direction or significance of
any change in group inference.

Next steps are a visual check of the corrected case and explicit propagation of
the corrected sm-5 input through affected UGR activation/PPI L2/L3 models and
extractions. Existing paper/group results have not yet been recomputed. Event
scanner-zero and original acquisition provenance remain investigation items.

For the insula analysis, use the completed unsmoothed sources and explicitly
replace this run's old sm-5 residual input for the smoothed comparison. Update
ICA input preparation to retain per-block acquisition TR and handle mixed TRs
without silently relabeling/interpolating the corrected data. Then run left,
right and bilateral ICA at dimensions 3/5/7/automatic for both residual scales.
Keep normalization policy consistent within a comparison; do not combine the
successful native-VN pilot with unrelated old normalization fits to declare a
complete homogeneous analysis. Older ICA fits containing the uncorrected run
remain preliminary and must be tracked as superseded when replacements run.
