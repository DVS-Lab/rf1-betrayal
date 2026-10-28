# Resume the existing queue with 64 slots

The [linux2 progress record](../linux2-trust-progress-20261028T195956Z-3cc48c/README.md), committed as 937ccfca, reported 96/108 successful fMRIPrep runs, 96 FEAT starts, 88 residuals and 88 numerical QC reports. All 88 terminal failures were at the subsequent reference-selection check. The resource sample showed 128 logical CPUs, 251.5 GiB total RAM and 76.7 GiB available; CPU idle was 94% and I/O wait 0% in its current one-second interval. Eight FEAT process trees were still active. This is a dated snapshot, not a live completion claim.

## Reference evidence

The three sampled actual logs explicitly list the run-2 SBRefs selected by workflow construction. For example, sub-10402's [log lines 87–88](../linux2-trust-progress-20261028T195956Z-3cc48c/sub-10402/reference-log-excerpt.txt) name all four run-2 SBRefs. The report figure's before/after SBRef line is absent; the checker had depended on that line. This produced an unknown-verification status after otherwise completed FEAT/QC. It does not justify repeating those fits.

The corrected checker still prefers the original before/after report evidence when available. Its additional route requires both an explicit selection-log list and a verified input view containing exactly the expected BOLD files for one run. Every listed reference must match that same run and exist in the view. Additional BOLD runs, mismatched reference names or missing evidence remain unverified. No filter, data, nuisance regressor or model is changed.

## Handoff and resources

The user requested a higher queue limit before launching the earlier 16-slot proposal. `start_resume_trust_stock_filter.sh` now launches `resume_trust_stock_filter.py` with **64 worker slots** and at most two new fMRIPrep launches. Set `INSULA_JOBS=32` when invoking the launcher to choose 32 instead. Existing FEAT jobs count within the worker limit. New stages also reserve RAM against MemAvailable: 28 GiB for fMRIPrep or 4 GiB for FEAT/QC, with another 12 GiB kept available. These reservations can reduce simultaneous new launches; 64 is a ceiling, not a promise of 64 active compute jobs throughout. The dated snapshot had only 20 runs without residuals remaining, plus metadata verification of existing fits.

The replacement waits for a boundary where the original Python scheduler's direct live children are exclusively FEAT commands for this output root. If a container, preflight or other child is active, the scheduler is immediately continued and checked again later. Before inspection, every scheduler thread must be stopped. The original scheduler alone receives termination; its FEAT children keep running and logging directly to their files. No fMRIPrep or FEAT compute process is terminated. The new scheduler records process start ticks, tracks inherited FEAT jobs and checks for other writers in each case before acting. A per-output advisory lock prevents a second replacement scheduler.

Successful fMRIPrep is reused. Existing residuals and FEAT fits are reused. Completed numerical QC is copied after small-file model/confound checks, header checks, a complete prior comparison report, matching residual path and file modification-time checks. When numerical QC needs repeating, an existing residual fit is still retained. A partial FEAT without residuals and without an active writer is preserved and reported for inspection, rather than overwritten. Unfinished fMRIPrep may resume in its existing work directory; fresh anatomical processing remains forbidden by the original cache preflight/normalization guard.

New records save the previous state, handoff identities, reference evidence, anatomy preflight, reused QC, successful input registry and any errors. They are pushed automatically when finished. No ICA or exclusions run.

## Local validation

The relevant 53-test suite passed before the concurrency increase. The updated nine-test resume suite also passes, including validation of 32/64 slots, tests against the three actual archived selection logs, rejection of mixed-run input views and wrong-run references, scheduler-only signalling, waiting through an active container, PID identity checks, RAM reservations, reuse of existing numerical QC, and an end-to-end simulated resume that never invokes preprocessing or FEAT for a completed case. Static Python checks and Bash syntax checking pass. Linux `/proc` discovery and handoff still require execution on linux2; local handoff tests simulate the process table and signals. The runtime guards validate ownership, paths, identities, stopped threads and direct child commands before termination.

```bash
cd /ZPOOL/data/projects/rf1-betrayal
git pull --ff-only origin insula-ica
bash code/insula_ica/start_resume_trust_stock_filter.sh
```
