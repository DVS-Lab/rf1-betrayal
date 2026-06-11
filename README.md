# Betrayal: Ultimatum Game and Trust Game Data and Analyses

This repository contains analysis code and selected derivative outputs for the manuscript:

> **Linking dACC Responses to Unreciprocated Trust to Brain and Behavior in the Ultimatum Game**

The hypotheses and analysis plans were pre-registered on the [Open Science Framework](https://osf.io/z6m35/overview).

The purpose of this repository is to support transparency and reproducibility. With access to the required raw behavioral files, preprocessed neuroimaging derivatives, templates, masks, and covariate files, users should be able to reproduce the primary behavioral and fMRI analyses reported in the manuscript.

> **Repository status:** this is not a complete standalone data release. Raw imaging data, large intermediate files, and some locally stored inputs are not tracked in git. Several scripts currently contain project-specific paths and should be reviewed before running on a new system.

---

## Analysis overview

The repository supports analyses for two fMRI tasks:

- **Trust Game** (`trust`): model `01`
- **Ultimatum Game** (`ugr`): model `3`, the updated model using the revised trial timing / decision-phase specification

The main analysis workflow is:

1. Identify participants with usable Trust Game and Ultimatum Game behavioral and neuroimaging data.
2. Convert raw task logs to BIDS-style `events.tsv` files.
3. Convert BIDS events into FSL-compatible three-column EV files.
4. Run behavioral, EV-file, and MRIQC-based exclusions.
5. Generate the final subject list used for fMRI analyses.
6. Run first-level, second-level, and third-level FSL analyses.
7. Run behavioral analyses in the associated R Markdown file.

---

## Repository organization

Expected repository structure:

```text
.
├── README.md
├── bids/                         # Generated BIDS-style events files; not necessarily complete/tracked
├── code/                         # Analysis and pipeline scripts
├── derivatives/
│   └── fsl/                      # FSL EV files and selected FEAT outputs
├── masks/                        # Seed/ROI masks used for PPI or ROI analyses
└── templates/                    # FSL .fsf templates for L1/L2/L3 models
```

Important notes:

- The repository does **not** include a complete BIDS dataset.
- Large files and raw imaging data are intentionally not tracked in git.
- Some derivative files may be tracked for transparency, but larger outputs can be regenerated if the required inputs are available.
- Several scripts infer the project root from the location of the `code/` directory, but others still use hard-coded local paths.

---

## Computing environment and dependencies

Users are expected to be comfortable working in a Linux command-line environment.

Required software:

- **Linux / Unix-like shell environment**
- **Bash**
- **FSL**, including `feat`, `fslnvols`, and `fslmeants`
- **MATLAB**, for converting raw behavioral logs to BIDS events
- **R**, for behavioral QC and subject-list generation

Required R packages:

```r
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
```

Additional environment assumptions:

- `FSLDIR` must be defined.
- The helper command `zeropad` is used in several scripts and must be available on the system path.
- The fMRI scripts assume fMRIPrep 24 derivatives and TedanaPlusConfounds files are available outside this repository.
- FSL `.fsf` templates must be present in `templates/`.
- Seed masks for PPI analyses must be present in `masks/` using the naming convention `seed-<seedname>.nii.gz`.

---

## Required inputs not stored in git

To reproduce the analyses, users need access to the following files or directories.

### Preprocessed fMRI data

The Level 1 scripts expect fMRIPrep outputs in the following form:

```text
/ZPOOL/data/projects/rf1-sra-linux2/derivatives/fmriprep-24/sub-<sub>/ses-01/func/
  sub-<sub>_ses-01_task-<trust|ugr>_run-<run>_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
```

### Confound regressors

The Level 1 scripts expect TedanaPlusConfounds files in the following form:

```text
/ZPOOL/data/projects/rf1-sra-linux2/derivatives/fsl/confounds_tedana-24/sub-<sub>/
  sub-<sub>_ses-01_task-<trust|ugr>_run-<run>_desc-TedanaPlusConfounds.tsv
```

If these confounds are missing, the Level 1 scripts exit before running FEAT.

### Raw behavioral logs

The conversion scripts expect raw task logs from the original RF1 project directories.

Trust Game logs:

```text
/ZPOOL/data/projects/rf1-sra/stimuli/Scan-Investment_Game/logs/
```

Ultimatum Game logs:

```text
/ZPOOL/data/projects/rf1-sra/stimuli/Scan-Lets_Make_A_Deal/logs/
```

### Subject-level covariate and QC files

`get_master_sublist.R` requires:

```text
Age.csv
AQ Data.csv
mriqc-metrics_allTasks_n299_ses-01.csv
behavioral_qc_dualtask.csv
ev_qc_dualtask.csv
subject_list_n225.txt
```

The script writes the final analysis subject list to:

```text
sublist_n132.txt
```

---

## Script inventory

### Subject-list and QC scripts

| Script | Purpose | Main output |
|---|---|---|
| `get_n225_sublist.sh` | Finds participants with both Trust and UGR preprocessed BOLD files and behavioral logs, then keeps the first 225 valid subjects. | `subject_list_all.txt`, `subject_list_n225.txt`, `exclusion_log.tsv`, `qc_table.tsv` |
| `check_EVfiles.sh` | Checks whether expected Trust and UGR EV files exist and are non-empty for each subject/run. | `ev_qc_dualtask.csv` |
| `check_dualtask_validtrials.R` | Applies the preregistered behavioral missing-trial exclusion. Trust requires at least 32/42 valid trials; UGR requires at least 36/48 valid trials. | `behavioral_qc_dualtask.csv` |
| `get_master_sublist.R` | Combines behavioral QC, EV QC, age exclusion, AQ availability, and MRIQC run-level exclusions into a final run-aligned analysis sample. | `master_subject_table_full.csv`, `sublist_n132.txt` |

### Event-conversion scripts

| Script | Purpose | Main output |
|---|---|---|
| `convertTrust_BIDS.m` | Converts raw Trust Game logs to BIDS-style events files. | `bids/sub-*/func/sub-*_task-trust_run-*_events.tsv` |
| `run_BIDSto3colTRUST.sh` | Loops over Trust Game BIDS event files and calls `BIDSto3colTRUST.sh` to generate FSL EV files. | `derivatives/fsl/EVfiles/sub-*/trust/` |
| `convertUGR_BIDS.m` | Converts raw Ultimatum Game logs to BIDS-style events files. | `bids/sub-*/func/sub-*_task-ugr_run-*_events.tsv` |
| `run_gen3colfilesUGR.sh` | Loops over subjects and calls `gen3colfilesUGR.sh` to generate UGR FSL EV files. | `derivatives/fsl/EVfiles/sub-*/ugr/model-3/` |

### fMRI analysis scripts

| Script | Purpose | Main output |
|---|---|---|
| `L1stats-trust.sh` | Runs Trust Game subject/run-level FEAT models. Supports activation and seed-based PPI modes. | `derivatives/fsl/sub-*/ses-01/*.feat` |
| `L1stats-ugr.sh` | Runs UGR subject/run-level FEAT models. Supports activation and seed-based PPI modes. | `derivatives/fsl/sub-*/ses-01/*.feat` |
| `run_L1stats-trust.sh` | Wrapper for Trust Game Level 1 models. | Multiple `.feat` directories |
| `run_L1stats-ugr.sh` | Wrapper for UGR Level 1 models. Currently configured for the `pTPJ` PPI seed. | Multiple `.feat` directories |
| `L2stats-trust.sh` | Combines Trust Game runs within subject. | `derivatives/fsl/sub-*/ses-01/*.gfeat` |
| `L2stats-ugr.sh` | Combines UGR runs within subject. | `derivatives/fsl/sub-*/ses-01/*.gfeat` |
| `run_L2stats-trust.sh` | Wrapper for Trust Game Level 2 models. | Multiple `.gfeat` directories |
| `run_L2stats-ugr.sh` | Wrapper for UGR Level 2 models. Currently configured for `ppi_seed-pTPJ`. | Multiple `.gfeat` directories |
| `L3stats-trust.sh` | Runs Trust Game group-level FEAT analyses. | Group-level `.gfeat` directories |
| `L3stats-ugr.sh` | Runs UGR group-level FEAT analyses. | Group-level `.gfeat` directories |
| `run_L3stats-trust.sh` | Wrapper for selected Trust Game Level 3 contrasts. Currently configured for cope 10, `rec-def`. | Group-level Trust output |

---

## Running the analyses

The scripts are intended to be run from the `code/` directory unless otherwise noted.

```bash
cd code
```

### 1. Generate the initial dual-task subject list

```bash
bash get_n225_sublist.sh
```

This creates:

```text
subject_list_all.txt
subject_list_n225.txt
exclusion_log.tsv
qc_table.tsv
```

`subject_list_n225.txt` is used as an input to later QC and event-conversion steps.

### 2. Convert raw Trust Game logs to BIDS events and FSL EV files

```bash
matlab -batch "run('convertTrust_BIDS.m')"
bash run_BIDSto3colTRUST.sh
```

`convertTrust_BIDS.m` writes BIDS-style `events.tsv` files. `run_BIDSto3colTRUST.sh` then converts those event files into FSL-compatible three-column EV files.

### 3. Convert raw UGR logs to BIDS events and FSL EV files

```bash
matlab -batch "run('convertUGR_BIDS.m')"
bash run_gen3colfilesUGR.sh
```

`convertUGR_BIDS.m` writes BIDS-style `events.tsv` files. `run_gen3colfilesUGR.sh` then converts those event files into FSL-compatible EV files for UGR model 3.

### 4. Check EV-file completeness

```bash
bash check_EVfiles.sh
```

This writes:

```text
ev_qc_dualtask.csv
```

The EV QC file records whether expected Trust and UGR run-level EV files are present and non-empty.

### 5. Apply behavioral valid-trial exclusions

```bash
Rscript check_dualtask_validtrials.R
```

This writes:

```text
behavioral_qc_dualtask.csv
```

The behavioral QC script applies the preregistered threshold requiring no more than 25% missing trials:

- Trust Game: at least 32 valid trials out of 42
- UGR: at least 36 valid trials out of 48

A subject is marked task-valid if at least one run passes the threshold for that task.

### 6. Generate the final analysis subject list

```bash
Rscript get_master_sublist.R
```

This script combines:

- the initial `subject_list_n225.txt` list,
- behavioral valid-trial QC,
- EV-file QC,
- age exclusion (`sub_age < 55`),
- AQ data availability,
- run-level MRIQC exclusions based on FD and tSNR outlier thresholds,
- run-aligned task validity.

It writes:

```text
master_subject_table_full.csv
sublist_n132.txt
```

If the subject list has Windows line endings, clean it before using it in Bash wrappers:

```bash
sed -i 's/\r$//' sublist_n132.txt
```

### 7. Run Level 1 fMRI models

Before running these wrappers, confirm that each wrapper uses the intended subject list and analysis type. Some wrappers currently point to older subject lists or specific PPI models.

Trust Game:

```bash
bash run_L1stats-trust.sh
```

UGR:

```bash
bash run_L1stats-ugr.sh
```

The Level 1 scripts:

- read fMRIPrep preprocessed BOLD files,
- read TedanaPlusConfounds files,
- generate subject/run-specific `.fsf` files from templates,
- run FSL FEAT,
- apply identity registration matrices for pre-normalized fMRIPrep data,
- delete selected large intermediate files after successful model estimation.

### 8. Run Level 2 fMRI models

Trust Game:

```bash
bash run_L2stats-trust.sh
```

UGR:

```bash
bash run_L2stats-ugr.sh
```

The Level 2 scripts combine run-level outputs within subject. Activation models expect 17 copes. PPI models add one contrast for the physiological regressor, producing 18 copes.

### 9. Run Level 3 fMRI models

Trust Game:

```bash
bash run_L3stats-trust.sh
```

UGR group-level models can be run by calling `L3stats-ugr.sh` directly or by adding a corresponding `run_L3stats-ugr.sh` wrapper.

Example direct call:

```bash
bash L3stats-ugr.sh <cope_number> <cope_name> <analysis_type>
```

The Level 3 scripts are currently configured for:

- final sample size: `N=132`
- covariate model: `full`
- FSL group model: `flame1`

### 10. Run behavioral analyses

Behavioral analyses are run from the associated R Markdown file:

```text
UG-Trust Behavioral Analyses.Rmd
```

Open or render this file in RStudio or from the command line after verifying that the input paths match the local repository structure.

---

## Reproducibility notes

This repository is intended to reproduce the reported analyses when paired with the required external data. It is not yet fully portable. In particular:

- Several scripts use hard-coded paths from the original analysis environment.
- Some wrapper scripts are configured for specific analysis types, seeds, contrasts, or subject lists.
- Helper scripts called by the wrappers must be present in `code/`.
- The `.fsf` templates in `templates/` must match the model names, contrast numbers, and placeholder names expected by the Bash scripts.
- The final subject list is generated by combining behavioral QC, EV QC, MRIQC, age, AQ availability, and run-alignment checks.

Before using this repository outside the original analysis environment, review all paths and wrapper settings carefully.

---


## Contact

Questions about this repository or the associated analyses can be directed to the manuscript authors.
