# Standardized task-bar review

All 490 runs and 2940 small OLS fits completed, generating eight original figures
and 150 plotted group bar estimates. Execution and Git publication reported no
errors. This does not establish numerical adequacy of every fitted coefficient.

## Trust: bilateral IC3 and IC4

The readable Trust-only figures in `review/` preserve every participant (N=132),
all run estimates and the original pointwise 95% CIs. No exclusion, clipping,
robust estimator or re-fit was used to generate them. They isolate the Trust
panel because one UGR estimate expands the common axis by orders of magnitude.

Joint-BOLD ICA basis, standardized reciprocate-minus-defect contrasts:

| Partner | IC3 mean [95% CI] | IC4 mean [95% CI] |
| --- | --- | --- |
| C | 0.01055 [-0.00173, 0.02283] | -0.01056 [-0.01933, -0.00180] |
| F | 0.02834 [0.01661, 0.04008] | -0.01498 [-0.02416, -0.00581] |
| S | 0.01231 [0.00054, 0.02409] | -0.00986 [-0.01917, -0.00055] |

IC3 means favor reciprocity, most clearly for F; its C interval spans zero. IC4
means favor defection for all three partners. Both ICA discovery bases show
closely similar directions and magnitudes. This is a same-cohort sensitivity
comparison, not independent replication. CIs are pointwise, not adjusted for the
150 plotted estimates. A positive-versus-negative pattern is descriptive and
is not a formal IC3-by-condition interaction test or evidence of unique encoding.

## UGR: numerical/design review required

Exactly one run accounts for every run-level coefficient/contrast with absolute
magnitude >5: `sub-11031_ses-01_task-ugr_run-1`. The affected estimates are
`social_high_pmod` and the high/low-averaged `ugr_social`, across all six map
bases and all five components (60 flagged rows). For bilateral joint-BOLD IC3,
the run's social-high standardized beta is 545893.924; for IC4 it is -99297.153.
These become group social means 1033.895 and -188.058 after within-participant
and group averaging. Their sizes reflect that one fit, not a credible basis for
interpreting a large population task response. Source matrices must be audited
before deciding the cause or an appropriate repair. No participant or run has
been removed. UGR interpretation is held pending that audit.

The saved original QC reports complete design preservation, exact original
nuisance identity and 100% insula coverage for this run. That does not rule out
collinearity, event coding problems or inadequate numerical rank handling.
`audit_network_task_design.py --publish` reads cached standardized matrices and
small task-event/design files only. It reports conditioning, each task EV's
conditional residual energy, event-amplitude summaries, and diagnostic cutoff
sensitivity for both runs of sub-11031. It neither reads 4D BOLD nor changes
production estimates, exclusion rules or the original first-level analysis.
