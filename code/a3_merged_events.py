import os
import numpy as np
import pandas as pd
import logging

# Define directories
data_dir = '/Users/tur61139/Documents/GitHub/rf1-norms/derivatives/fsl/EVfiles'
subjects = [subj for subj in os.listdir(data_dir) if subj.startswith('sub-')]
log_file = os.path.join('/Users/tur61139/Documents/GitHub/rf1-norms/logs', "a3_mergeLog.txt")
logging.basicConfig(filename=log_file, level=logging.WARNING,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Runs and trial types to process
runs = ['run-1', 'run-2']
trial_types = ['nonsocial_high', 'nonsocial_low', 'social_high', 'social_low']

# TSV input data directory
tsv_data_dir = '/Users/tur61139/Documents/GitHub/rf1-norms/bids'

# ORIGINAL SCRIPT - MODIFIED TO INCLUDE MISSED TRIAL PROCESSING
for subject in subjects:
    logging.info(f"Processing subject: {subject}")
    subj_dir = os.path.join(data_dir, subject, 'ugr')
    output_dir = os.path.join(subj_dir, 'merged_events')
    os.makedirs(output_dir, exist_ok=True)

    for run in runs:
        # Check for missed trial input file
        missed_trial_file = os.path.join(subj_dir, f"{run}_missed_trial_constant.txt")
        if os.path.isfile(missed_trial_file):
            missed_trial_data = pd.read_csv(missed_trial_file, delim_whitespace=True, header=None)
            if missed_trial_data.shape[1] < 3:  # Corrected condition
                logging.warning(f"Malformed missed trial input file: {missed_trial_file}, skipping.")
            else:
                missed_onsets = missed_trial_data.iloc[:, 0] - 2
                missed_durations = missed_trial_data.iloc[:, 1]
                missed_values = np.ones(len(missed_onsets))
                missed_output = np.column_stack([missed_onsets, missed_durations, missed_values])
                missed_output_file = os.path.join(output_dir, f"{run}_missed_trial.txt")
                np.savetxt(missed_output_file, missed_output, fmt="%.6f\t%.6f\t%.1f", delimiter="\t")
                print(f"  Missed trial file created for {run}.")

        for trial_type in trial_types:
            dec_file = os.path.join(subj_dir, f"{run}_dec_{trial_type}_constant.txt")
            if not os.path.isfile(dec_file):
                logging.warning(f"  Missing file: {dec_file}, skipping {trial_type} for {run}.")
                continue

            dec_data = pd.read_csv(dec_file, delim_whitespace=True, header=None)
            if dec_data.shape[1] < 3:  # Corrected condition
                logging.warning(f"  Malformed decision input file: {dec_file}, skipping.")
                continue

            onset_times = dec_data.iloc[:, 0] - 2
            durations = dec_data.iloc[:, 1] + 2
            constant_values = np.ones_like(dec_data.iloc[:, 2])
            demeaned_values = dec_data.iloc[:, 2] - dec_data.iloc[:, 2].mean()

            constant_output = np.column_stack([onset_times, durations, constant_values])
            constant_file = os.path.join(output_dir, f"{run}_{trial_type}_constant.txt")
            np.savetxt(constant_file, constant_output, fmt="%.6f\t%.6f\t%.1f", delimiter="\t")

            pmod_output = np.column_stack([onset_times, durations, demeaned_values])
            pmod_file = os.path.join(output_dir, f"{run}_{trial_type}_pmod.txt")
            np.savetxt(pmod_file, pmod_output, fmt="%.6f", delimiter="\t")

            print(f"  Processed {trial_type} for {run}.")

# AMENDMENTS - ONLY ADDS.TSV-BASED UPDATES TO OUTPUT FILES
for subject in subjects:
    print(f"Amending files for subject: {subject}")
    subj_dir = os.path.join(data_dir, subject, 'ugr')
    output_dir = os.path.join(subj_dir, 'merged_events')

    for run in runs:
        tsv_file = os.path.join(tsv_data_dir, subject, 'func', f"{subject}_task-ugr_{run}_events.tsv")
        if not os.path.isfile(tsv_file):
            logging.warning(f"Missing TSV file: {tsv_file}")
            continue

        # Load the.tsv file
        tsv_data = pd.read_csv(tsv_file, sep="\t", header=None)

        # Ensure numeric conversion for relevant columns
        tsv_data.iloc[:, 0] = pd.to_numeric(tsv_data.iloc[:, 0], errors='coerce')  # Onsets
        tsv_data.iloc[:, 1] = pd.to_numeric(tsv_data.iloc[:, 1], errors='coerce')  # Durations
        tsv_data.iloc[:, 2] = pd.to_numeric(tsv_data.iloc[:, 2], errors='coerce')  # Decision values


        # Group cues with the same onset time
        grouped_cues = {}
        for i, row in tsv_data.iterrows():
            if "cue" in row:  # Check for "cue" in the correct column
                onset_time = row  # Extract onset time from the correct column
                if onset_time not in grouped_cues:
                    grouped_cues[onset_time] = []
                grouped_cues[onset_time].append(row)

        # Pair grouped cues with dec/missed_trial
        paired_data = []
        cue_indices = list(grouped_cues.keys())
        for i, (onset_time, cue_rows) in enumerate(grouped_cues.items()):
            next_row = None
            if i + 1 < len(cue_indices):
                next_onset = cue_indices[i + 1]
                next_row = tsv_data[tsv_data.iloc[:, 0] == next_onset].iloc[0] if not tsv_data[tsv_data.iloc[:, 0] == next_onset].empty else None
            if next_row is not None and ("dec" in next_row or "missed_trial" in next_row):
                for cue_row in cue_rows:
                    cue_onset = cue_row  # Corrected
                    event_onset = next_row  # Corrected
                    event_duration = next_row  # Corrected
                    if pd.notna(cue_onset) and pd.notna(event_onset) and pd.notna(event_duration):
                        duration = (event_onset + event_duration) - cue_onset
                        value = 1.0 if "missed_trial" in next_row else next_row.iloc[2] - tsv_data.iloc[:, 2].mean()
                        paired_data.append([cue_onset, duration, value])

        # Save paired data to relevant files
        for trial_type in trial_types + ["missed_trial"]:  # Include missed_trial
            constant_file = os.path.join(output_dir, f"{run}_{trial_type}_constant.txt")

            if os.path.isfile(constant_file):
                constant_data = [row for row in paired_data if trial_type in row]
                constant_data = np.array(constant_data)

                if constant_data.size > 0:
                    if trial_type == "missed_trial":
                        # Add a column of 1.0s for missed trials
                        ones_col = np.ones((constant_data.shape[0], 1))
                        constant_data = np.hstack((constant_data, ones_col))
                        fmt_str = "%.6f\t%.6f\t%.1f"  # 3 columns for missed trials
                    else:
                        fmt_str = "%.6f\t%.6f\t%.1f"  # 3 columns for other trials

                    np.savetxt(constant_file, constant_data, fmt=fmt_str, delimiter="\t")
                    print(f"  Amended {constant_file}.")
                else:
                    logging.warning(f"No data to write to {constant_file} for subject: {subject}, run: {run}, trial_type: {trial_type}")

            pmod_file = os.path.join(output_dir, f"{run}_{trial_type}_pmod.txt")
            if os.path.isfile(pmod_file) and trial_type!= "missed_trial":  # Exclude missed_trial
                pmod_data = [row for row in paired_data if f"cue_{trial_type}" in row]
                pmod_data = np.array(pmod_data)

                if pmod_data.size > 0:
                    np.savetxt(pmod_file, pmod_data, fmt="%.6f", delimiter="\t")
                    print(f"  Amended {pmod_file}.")
                else:
                    logging.warning(f"No data to write to {pmod_file} for subject: {subject}, run: {run}, trial_type: {trial_type}")

print("Processing complete.")