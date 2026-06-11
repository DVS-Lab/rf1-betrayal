library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# ===============================
# CONFIGURATION
# ===============================

base_dir <- "D:/Postdoc/UG-Trust/code"

master_list_file  <- file.path(base_dir, "subject_list_n225.txt")
behavioral_file   <- file.path(base_dir, "behavioral_qc_dualtask.csv")
ev_qc_file        <- file.path(base_dir, "ev_qc_dualtask.csv")
age_file          <- file.path(base_dir, "Age.csv")
aq_file           <- file.path(base_dir, "AQ Data.csv")
mriqc_file        <- file.path(base_dir, "mriqc-metrics_allTasks_n299_ses-01.csv")

# ===============================
# LOAD MASTER LIST
# ===============================

master_ids <- read_lines(master_list_file) %>%
  str_remove("^sub-") %>%
  as.numeric()

master_df <- data.frame(sub_id = master_ids)

cat("Starting N:", nrow(master_df), "\n")

# ===============================
# MERGE BEHAVIORAL QC (RENAMED)
# ===============================

behavioral_qc <- read_csv(behavioral_file, show_col_types = FALSE) %>%
  rename(
    trust_beh_run1_ok = trust_run1_ok,
    trust_beh_run2_ok = trust_run2_ok,
    trust_beh_valid   = trust_valid,
    ugr_beh_run1_ok   = ugr_run1_ok,
    ugr_beh_run2_ok   = ugr_run2_ok,
    ugr_beh_valid     = ugr_valid
  )

master_df <- master_df %>%
  left_join(behavioral_qc, by = "sub_id")

# ===============================
# MERGE EV QC (RENAMED)
# ===============================

ev_qc <- read_csv(ev_qc_file, show_col_types = FALSE) %>%
  rename(
    trust_ev_run1_ok = trust_run1_ok,
    trust_ev_run2_ok = trust_run2_ok,
    trust_ev_status  = trust_valid_runs,
    ugr_ev_run1_ok   = ugr_run1_ok,
    ugr_ev_run2_ok   = ugr_run2_ok,
    ugr_ev_status    = ugr_valid_runs
  )

master_df <- master_df %>%
  left_join(ev_qc, by = "sub_id")

# ===============================
# APPLY AGE EXCLUSION (<55)
# ===============================

Age <- read_csv(age_file, show_col_types = FALSE)

Age <- Age %>%
  mutate(
    sub_id = as.numeric(str_extract(sub_id, "\\d+"))
  ) %>%
  filter(!is.na(sub_id)) %>%
  filter(redcap_event_name == "subject_informatio_arm_1") %>%
  select(sub_id, sub_age)

master_df <- master_df %>%
  left_join(Age, by = "sub_id") %>%
  filter(sub_age < 55) %>%
  distinct(sub_id, .keep_all = TRUE)

cat("After age exclusion:", nrow(master_df), "\n")

# ===============================
# APPLY AQ EXCLUSION IF MISSING DATA
# ===============================

AQ <- read_csv(aq_file, show_col_types = FALSE)

AQ <- AQ %>%
  mutate(
    sub_id = as.numeric(str_extract(sub_id, "\\d+"))
  ) %>%
  filter(!is.na(sub_id)) %>%
  filter(redcap_event_name == "subject_informatio_arm_1") %>%
  select(sub_id, score_aq_total)

master_df <- master_df %>%
  left_join(AQ, by = "sub_id") %>%
  filter(!is.na(score_aq_total)) %>%
  distinct(sub_id, .keep_all = TRUE)

cat("After AQ exclusion:", nrow(master_df), "\n")


# ===============================
# MRIQC PROCESSING
# ===============================

mriqc <- read_csv(mriqc_file, show_col_types = FALSE) %>%
  mutate(sub = as.numeric(sub)) %>%
  filter(sub %in% master_df$sub_id) %>%
  filter(task %in% c("ugr", "trust"))

mriqc_flagged <- mriqc %>%
  group_by(task) %>%
  mutate(
    fd_threshold =
      quantile(fd_mean, 0.75, na.rm = TRUE) +
      1.5 * IQR(fd_mean, na.rm = TRUE),
    
    tsnr_threshold =
      quantile(tsnr, 0.25, na.rm = TRUE) -
      1.5 * IQR(tsnr, na.rm = TRUE),
    
    remove_run =
      fd_mean > fd_threshold |
      tsnr < tsnr_threshold
  ) %>%
  ungroup()

valid_runs_df <- mriqc_flagged %>%
  filter(!remove_run)

task_summary <- valid_runs_df %>%
  group_by(sub, task) %>%
  summarise(
    remaining_runs = paste(sort(run), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(
    valid_runs = case_when(
      remaining_runs == "1,2" ~ "both",
      remaining_runs == "1"   ~ "run1",
      remaining_runs == "2"   ~ "run2"
    )
  )

mriqc_wide <- task_summary %>%
  select(sub, task, valid_runs) %>%
  pivot_wider(
    names_from = task,
    values_from = valid_runs,
    names_prefix = "mriqc_"
  ) %>%
  rename(sub_id = sub)

master_df <- master_df %>%
  left_join(mriqc_wide, by = "sub_id") %>%
  mutate(
    mriqc_trust = ifelse(is.na(mriqc_trust), "none", mriqc_trust),
    mriqc_ugr   = ifelse(is.na(mriqc_ugr), "none", mriqc_ugr)
  )

# ===============================
# EXPAND MRIQC INTO RUN FLAGS
# ===============================

master_df <- master_df %>%
  mutate(
    trust_mri_run1_ok = ifelse(mriqc_trust %in% c("both","run1"), 1, 0),
    trust_mri_run2_ok = ifelse(mriqc_trust %in% c("both","run2"), 1, 0),
    ugr_mri_run1_ok   = ifelse(mriqc_ugr %in% c("both","run1"), 1, 0),
    ugr_mri_run2_ok   = ifelse(mriqc_ugr %in% c("both","run2"), 1, 0)
  )

# ===============================
# FULL RUN VALIDITY (ALL QC TYPES)
# ===============================

master_df <- master_df %>%
  mutate(
    trust_run1_full_ok =
      trust_beh_run1_ok == 1 &
      trust_ev_run1_ok  == 1 &
      trust_mri_run1_ok == 1,
    
    trust_run2_full_ok =
      trust_beh_run2_ok == 1 &
      trust_ev_run2_ok  == 1 &
      trust_mri_run2_ok == 1,
    
    ugr_run1_full_ok =
      ugr_beh_run1_ok == 1 &
      ugr_ev_run1_ok  == 1 &
      ugr_mri_run1_ok == 1,
    
    ugr_run2_full_ok =
      ugr_beh_run2_ok == 1 &
      ugr_ev_run2_ok  == 1 &
      ugr_mri_run2_ok == 1
  )

# ===============================
# TASK-LEVEL VALIDITY (≥1 FULL RUN)
# ===============================

master_df <- master_df %>%
  mutate(
    trust_task_valid = trust_run1_full_ok | trust_run2_full_ok,
    ugr_task_valid   = ugr_run1_full_ok   | ugr_run2_full_ok
  )

# ===============================
# FINAL INCLUSION (RUN-ALIGNED)
# ===============================

master_df <- master_df %>%
  mutate(
    include_final = trust_task_valid & ugr_task_valid
  )

master_inclusion_df <- master_df %>%
  filter(include_final)

cat("Final N after ALL exclusions (run-aligned):",
    nrow(master_inclusion_df), "\n")


# ===============================
# EXPORT
# ===============================

write.table(
  master_inclusion_df$sub_id,
  file = "sublist_n132.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.csv(
  master_df,
  file = "master_subject_table_full.csv",
  row.names = FALSE
)

cat("Master table written.\n")
