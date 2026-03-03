# Betrayal data: Analysis and Statistical Modeling
This repository contains code and selected derivative outputs for analyses related to our manuscript, **"Linking dACC Responses to Unreciprocated Trust to Brain and Behavior in the Ultimatum Game."** It is designed as an analysis-focused companion repository (primarily FSL modeling and ROI extraction) rather than a full raw-data archive. Analysis plans were pre-registered on [Open Science Framework][osf].


## A few prerequisites and recommendations
- Understand BIDS-style derivatives and be comfortable navigating Linux
- Install [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
- Use a shell environment with `bash`
- Ensure `feat`, `fslmaths`, `fslmeants`, `fslnvols`, and `zeropad` are available in your PATH
- Ideally have access to both Linux2 and HPC storage if you are reproducing the original compute workflow
- Expect to edit hard-coded absolute paths before running scripts on a new system


## Notes on repository organization and files
- This repository does **not** contain a complete raw BIDS dataset.
- Large imaging files are primarily expected in external storage locations and may not be tracked in git.
- Tracked folders and their contents:
  - `code`: run scripts, level-specific FSL scripts, extraction scripts, subject lists, and QC helper files
  - `templates`: FSL design templates (`.fsf`) and related design artifacts (`.con`, `.mat`, `.grp`, `.png`, `.ppm`)
  - `masks`: ROI and intersection masks used for PPI and signal extraction
  - `derivatives`: selected outputs (including extracted ROI signals) for transparency and reproducibility


## Downloading Data and Running Analyses
```
# Clone repository
git clone https://github.com/DVS-Lab/rf1-betrayal.git
cd rf1-betrayal
```

### Step 1: Confirm and edit path variables before running anything
Most analysis scripts assume data live in project-specific absolute paths (for example `/ZPOOL/data/projects/rf1-sra-linux2` and `/ZPOOL/data/projects/rf1-betrayal`).

At minimum, review and update the following scripts for your environment:
```
code/L1stats-ugr.sh
code/L1stats-trust.sh
code/L3stats-ugr.sh
code/L3stats-trust.sh
code/get_n225_sublist.sh
code/getfilepaths-trust.sh
code/transfer_EVfiles.sh
code/extract_trust-dACC.sh
```

### Step 2: Confirm required inputs exist
Before running models, verify availability of:
- Preprocessed task BOLD data (MNI space)
- Tedana/confound TSV files used by L1 models
- EV files in the expected derivatives EV directory
- Subject lists (for example `sublist_DD128.txt`) aligned with available data

Typical external dependencies referenced by scripts include:
```
/ZPOOL/data/projects/rf1-sra-linux2/derivatives/fmriprep-24
/ZPOOL/data/projects/rf1-sra-linux2/derivatives/fsl/confounds_tedana-24
/ZPOOL/data/projects/rf1-betrayal/derivatives/fsl/EVfiles
/ZPOOL/data/projects/rf1-sra/stimuli
```

### Step 3: Run UGR Level 1 analyses
```
cd code
bash run_L1stats-ugr.sh
```
This wrapper loops through `sublist_DD128.txt`, both runs, and both analysis modes currently configured in-script (`act` and `ppi` with `aIns`).

### Step 4: Run UGR Level 2 analyses
```
cd code
bash run_L2stats-ugr.sh
```
This combines run-level outputs into subject-level fixed-effects results.

### Step 5: Run UGR Level 3 analyses
```
cd code
bash run_L3stats-ugr.sh
```
This loops across predefined cope sets and runs group-level FEAT analyses using configured L3 templates.

### Step 6: Run Trust Level 1 analyses
```
cd code
bash run_L1stats-trust.sh
```
Current defaults are configured for activation analysis (`ppi=0`).

### Step 7: Run Trust Level 2 and Level 3 analyses
```
cd code
bash run_L2stats-trust.sh
bash run_L3stats-trust.sh
```
These scripts produce subject-level and group-level trust outputs using the currently configured model/type/template settings.

### Step 8: Run ROI extraction workflows (optional)
If you need ROI-level extracted signals from task contrasts, use scripts such as:
```
cd code
bash extract_ugr_AIns.sh
bash extract_ugr_AIns_social.sh
bash extract_trust_AIns.sh
bash extract_trust_AIns_social.sh
```
Additional extraction scripts exist for specialized/legacy workflows (`extract_trust-dACC.sh`, `extract_fevsxpanas-n-change.sh`) and may require path edits.

### Step 9: Generate or refresh cohort/QC helper outputs (optional)
To regenerate dual-task inclusion and QC tables:
```
cd code
bash get_n225_sublist.sh
```
Outputs include:
- `subject_list_all.txt`
- `subject_list_n225.txt`
- `qc_table.tsv`
- `exclusion_log.tsv`


## Additional implementation notes
- Wrapper scripts use background execution (`&`) and a process-count throttle (`NCORES`) to manage concurrency.
- Because jobs are submitted in the background, verify completion of one stage before launching the next stage.
- Several scripts remove large intermediate FEAT files after completion to conserve disk space.
- L3 scripts assume specific covariate names, model numbers, and template files; verify these before reruns.


## Acknowledgments
This work was supported, in part, by grants from the National Institutes of Health.

[osf]: https://osf.io/z6m35/overview
