import pandas as pd
import os

# Paths
sublist_path = "/Users/tur61139/Documents/GitHub/rf1-norms/code/sublist_validbehUGR.txt"
log_path = "/Users/tur61139/Documents/GitHub/rf1-norms/logs/a4_mergeLog.txt"
output_dir_template = "/Users/tur61139/Documents/GitHub/rf1-norms/derivatives/fsl/EVfiles/sub-{subject}/ugr/model-3"
raw_data_template = "/Users/tur61139/Documents/GitHub/rf1-sra/stimuli/Scan-Lets_Make_A_Deal/logs/{subject}/sub-{subject}_task-ultimatum_run-{run}_raw.csv"

# Read subject list
with open(sublist_path, "r") as f:
    subjects = [line.strip() for line in f]

# Initialize log
log = []

def process_subject(subject, run):
    if subject == "10817" and run == "2":
        log.append(f"Skipping {subject} run-{run}: Excluded due to missing data.")
        return
    
    raw_csv_path = raw_data_template.format(subject=subject, run=str(int(run)-1))
    if not os.path.exists(raw_csv_path):
        log.append(f"Skipping {subject} run-{run}: No input file found.")
        return
    
    df = pd.read_csv(raw_csv_path)
    if "resp_onset" not in df.columns or "rt" not in df.columns or "Endowment" not in df.columns or "L_Option" not in df.columns or "R_Option" not in df.columns or "cue_Onset" not in df.columns or "decision_offset" not in df.columns:
        log.append(f"Skipping {subject} run-{run}: Missing necessary columns.")
        return
    
    df["rt"] = pd.to_numeric(df["rt"], errors="coerce")
    
    # Ensure 'decision_offset' and 'cue_Onset' are numeric
    if not pd.api.types.is_numeric_dtype(df["decision_offset"]):
        df["decision_offset"] = pd.to_numeric(df["decision_offset"], errors="coerce")
    if not pd.api.types.is_numeric_dtype(df["cue_Onset"]):
        df["cue_Onset"] = pd.to_numeric(df["cue_Onset"], errors="coerce")
    
    valid_rt = df.loc[df["rt"] != 999, "rt"]
    mean_rt = valid_rt.mean() if not valid_rt.empty else 0
    df["demeaned_rt"] = df["rt"] - mean_rt
    df["zero"] = 0
    df["duration"] = df["decision_offset"] - df["cue_Onset"]
    df["offer_amount"] = df["L_Option"] + df["R_Option"]

    def demean_offer(group):
        # Explicitly select non-grouping columns
        group["demeaned_offer"] = group["offer_amount"] - group["offer_amount"].mean()
        return group

    df = df.groupby(["Block", "Endowment"]).apply(demean_offer).reset_index(drop=True)

    output_dir = output_dir_template.format(subject=subject)
    os.makedirs(output_dir, exist_ok=True)
    
    for block in ["social", "nonsocial"]:
        for level in ["high", "low"]:
            event_file = f"run-{run}_{block}_{level}_constant.txt"
            pmod_file = f"run-{run}_{block}_{level}_pmod.txt"
            
            event_path = os.path.join(output_dir, event_file)
            pmod_path = os.path.join(output_dir, pmod_file)
            
            subset = df[(df["Block"] == (3 if block == "social" else 2)) & 
                        (df["Endowment"] == (32 if level == "high" else 16)) &
                        (df["rt"] != 999)]
            
            if not subset.empty:
                # Sort by 'cue_Onset' before writing
                subset = subset.sort_values(by='cue_Onset')
                subset[["cue_Onset", "duration"]].assign(one=1.0).to_csv(event_path, sep="\t", index=False, header=False, float_format="%.6f")
                subset[["cue_Onset", "duration", "demeaned_offer"]].to_csv(pmod_path, sep="\t", index=False, header=False, float_format="%.6f")
    
    missed_trials = df[df["rt"] == 999]
    if not missed_trials.empty:
        missed_trials_path = os.path.join(output_dir, f"run-{run}_missed_trial.txt")
        # Sort by 'cue_Onset' before writing
        missed_trials = missed_trials.sort_values(by='cue_Onset')
        missed_trials[["cue_Onset", "duration"]].assign(one=1.0).to_csv(missed_trials_path, sep="\t", index=False, header=False, float_format="%.6f")
    
    rt_constant_path = os.path.join(output_dir, f"run-{run}_rt_constant.txt")
    rt_pmod_path = os.path.join(output_dir, f"run-{run}_rt_pmod.txt")

    valid_df = df[df["rt"] != 999]
    
    if not valid_df.empty:
        # Sort by 'resp_onset' before writing
        valid_df = valid_df.sort_values(by='resp_onset')
        valid_df[["resp_onset"]].assign(zero=0, one=1.0).to_csv(rt_constant_path, sep="\t", index=False, header=False, float_format="%.6f")
        valid_df[["resp_onset", "zero", "demeaned_rt"]].to_csv(rt_pmod_path, sep="\t", index=False, header=False, float_format="%.6f")
    
    log.append(f"Processed {subject} run-{run}")

for subject in subjects:
    for run in ["1", "2"]:
        process_subject(subject, run)

with open(log_path, "w") as f:
    f.write("\n".join(log))
