# Redo functional inputs with the stock filter and reuse completed anatomy

The user directed us to redo the affected inputs and use the original fMRIPrep configuration. The run-list override in the generated SBRef filter was our error. The original wrapper and its filter are not edited.

## Source of authority

The snapshots in sources/ were retrieved from [rf1-sra-linux2/code/fmriprep.sh](https://github.com/DVS-Lab/rf1-sra-linux2/blob/main/code/fmriprep.sh) and [fmriprep_config.json](https://github.com/DVS-Lab/rf1-sra-linux2/blob/main/code/fmriprep_config.json) on 2026-10-28. The stock filter specifies datatype, suffix and permitted parts, with no run override. The rerun reads the actual linux2 filter, requires semantic agreement with this snapshot, and copies its exact bytes unchanged. It also records the original linux2 wrapper and its digest.

This remains the agreed **24.1.1 volumetric analysis** (T1w/MNI152NLin6Asym, no FreeSurfer reconstruction). It does not invoke the modern wrapper's default image/CIFTI/FreeSurfer workflow. Added orchestration is limited to selecting the required input files in a read-only BIDS view, routing outputs and work to fresh directories, allocating eight concurrent jobs, and supplying existing anatomical derivatives. No run list is injected into JSON queries. The stock default reference selection, including its BOLD-derived fallback when a corresponding SBRef is absent, is retained. The earlier fixed-seed additions are not part of this stock-default rerun; the anatomical results themselves are reused unchanged.

## Scope: 108 remaining Trust run-2 inputs

The completed repair audit identified 108 cross-run selections and two missing-log cases. The archived **sub-10478-v24-a** log resolves the latter: its run 1 used run 1's SBRef and its run 2 also used run 1's SBRef. This brings the original rebuilt set to 132 same-run and 109 cross-run references. Sub-10834's reviewed two-run repair is retained as an explicit exception. The remaining queue is **108 run-2 inputs from 108 participants**, processed eight at a time. Other run-1 inputs and all UGR inputs are not rerun.

## Avoid repeating anatomical ANTs registration

For each selected participant, stage only its completed **anatomical** derivatives as read-only links: native preprocessed T1w, brain mask, segmentation, three tissue-probability maps, and forward/reverse transforms for MNI152NLin6Asym and MNI152NLin2009cAsym. No old functional references, motion transforms or BOLD-to-T1w transforms are exposed to derivative reuse.

Before processing, the installed container must discover all required cached assets and construct an anatomical graph with no new normalization, skull-stripping, bias-correction, segmentation or T1w-template-estimation interfaces. Missing/ambiguous/incompatible assets stop that job. The actual run is also monitored for a request to prepare new anatomical normalization and is terminated if that occurs. Completion requires the actual log to confirm precomputed anatomical registrations were found.

This uses the documented [fMRIPrep 24.1.1 derivative-reuse mechanism](https://fmriprep.org/en/24.1.1/usage.html#bids-derivatives-reuse). Functional reference creation, distortion correction, registration, resampling, optimal combination and the original joint-model FEAT fit still need to run. Functional processing may use ANTs tools; the expensive T1w-to-template anatomical registration is the part being reused.

## Verification and records

The driver verifies the same 24.1.1 image digest, exact original nuisance-file values/hash, original full model matrices, grid/TR/time count, finite nonconstant residuals, actual run-matched reference selection (or the stock missing-SBRef fallback), and unchanged cached anatomical file identities. It saves a QC montage and left/right/bilateral coverage per run. New inputs remain separate until review; no exclusions or ICA run automatically. A failed-subject subset can be retried with `rerun_trust_references.py --subjects ... --publish` without redoing successful participants.

Run on linux2:

```bash
cd /ZPOOL/data/projects/rf1-betrayal
git pull --ff-only origin insula-ica
bash code/insula_ica/start_trust_stock_filter_rerun.sh
```

The launcher prints a short log command and runs in the background. Completion/failure records are committed and pushed automatically. Actual cached-anatomy reuse must still be established by the linux2 preflight/run; local tests validate the orchestration and stop conditions without executing MRI processing.

Local validation: all 39 tests across stock rerun, repair, rebuild, regeneration and BOLD controls pass. They include a real child-process termination check, unchanged filter bytes, selected-run input scoping, anatomy-only derivative exposure, and a simulated completed pipeline. Python static checks, Bash syntax checking and Git whitespace checking also pass. Container execution and MRI results remain to be verified on linux2.
