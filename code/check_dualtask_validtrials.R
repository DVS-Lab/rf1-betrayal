library(dplyr)
library(readr)
library(stringr)

# ===============================
# CONFIGURATION
# ===============================

base_dir <- "D:/Postdoc/UG-Trust/code"

trust_dir <- file.path(base_dir, "trust_logs")
ugr_dir   <- file.path(base_dir, "ugr_logs")

subject_list_file <- file.path(base_dir, "subject_list_n225.txt")

TRUST_TOTAL <- 42
TRUST_THRESHOLD <- ceiling(0.75 * TRUST_TOTAL)   # 32

UGR_TOTAL <- 48
UGR_THRESHOLD <- ceiling(0.75 * UGR_TOTAL)       # 36

# ===============================
# LOAD SUBJECT LIST
# ===============================

subjects <- read_lines(subject_list_file)
subjects <- str_remove(subjects, "^sub-")
subjects <- as.numeric(subjects)

# ===============================
# STORAGE
# ===============================

results <- list()

# ===============================
# LOOP THROUGH SUBJECTS
# ===============================

for (sub in subjects) {
  
  # -------------------------------
  # TRUST CHECK
  # -------------------------------
  
  trust_run1_ok <- 0
  trust_run2_ok <- 0
  
  for (r in 0:1) {
    
    file_path <- file.path(
      trust_dir,
      sub,
      sprintf("sub-%05d_task-trust_run-%d_raw.csv", sub, r)
    )
    
    if (!file.exists(file_path)) next
    
    T <- tryCatch(read_csv(file_path, show_col_types = FALSE),
                  error = function(e) NULL)
    
    if (is.null(T) || !"resp" %in% colnames(T)) next
    
    valid_trials <- sum(T$resp != 999, na.rm = TRUE)
    
    if (valid_trials >= TRUST_THRESHOLD) {
      if (r == 0) trust_run1_ok <- 1
      if (r == 1) trust_run2_ok <- 1
    }
  }
  
  trust_valid <- ifelse(trust_run1_ok == 1 | trust_run2_ok == 1, 1, 0)
  
  # -------------------------------
  # UGR CHECK
  # -------------------------------
  
  ugr_run1_ok <- 0
  ugr_run2_ok <- 0
  
  for (r in 0:1) {
    
    file_path <- file.path(
      ugr_dir,
      sub,
      sprintf("sub-%05d_task-ultimatum_run-%d_raw.csv", sub, r)
    )
    
    if (!file.exists(file_path)) next
    
    T <- tryCatch(read_csv(file_path, show_col_types = FALSE),
                  error = function(e) NULL)
    
    if (is.null(T) || !"resp" %in% colnames(T)) next
    
    valid_trials <- sum(T$resp %in% c(1,2), na.rm = TRUE)
    
    if (valid_trials >= UGR_THRESHOLD) {
      if (r == 0) ugr_run1_ok <- 1
      if (r == 1) ugr_run2_ok <- 1
    }
  }
  
  ugr_valid <- ifelse(ugr_run1_ok == 1 | ugr_run2_ok == 1, 1, 0)
  
  # -------------------------------
  # STORE
  # -------------------------------
  
  results[[as.character(sub)]] <- data.frame(
    sub_id = sub,
    trust_run1_ok = trust_run1_ok,
    trust_run2_ok = trust_run2_ok,
    trust_valid = trust_valid,
    ugr_run1_ok = ugr_run1_ok,
    ugr_run2_ok = ugr_run2_ok,
    ugr_valid = ugr_valid
  )
}

# ===============================
# COMBINE + WRITE
# ===============================

behavioral_qc <- bind_rows(results)

write_csv(behavioral_qc,
          file.path(base_dir, "behavioral_qc_dualtask.csv"))

cat("Behavioral QC complete.\n")
cat("Trust valid:", sum(behavioral_qc$trust_valid), "\n")
cat("UGR valid:", sum(behavioral_qc$ugr_valid), "\n")
