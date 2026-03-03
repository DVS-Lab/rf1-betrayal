# read subject list
subs <- scan("sublist_n132.txt", what = character())

# initialize output data frame
out <- data.frame(
  sub_id = subs,
  nonsocial_lvl1 = NA_real_,
  nonsocial_lvl2 = NA_real_,
  nonsocial_lvl3 = NA_real_,
  nonsocial_lvl4 = NA_real_,
  social_lvl1 = NA_real_,
  social_lvl2 = NA_real_,
  social_lvl3 = NA_real_,
  social_lvl4 = NA_real_,
  stringsAsFactors = FALSE
)

# helper function to read the correct file
read_meants <- function(folder, sub) {
  l2_file  <- file.path(folder, paste0("sub-", sub, "_L2_AIns_meants.txt"))
  run1_file <- file.path(folder, paste0("sub-", sub, "_run-1_AIns_meants.txt"))
  
  if (file.exists(l2_file)) {
    as.numeric(scan(l2_file, what = numeric(), quiet = TRUE))
  } else {
    as.numeric(scan(run1_file, what = numeric(), quiet = TRUE))
  }
}

# loop over subjects
for (i in seq_along(subs)) {
  s <- subs[i]
  
  out$nonsocial_lvl1[i] <- read_meants("AIns_cat1_ugr_meants", s)
  out$nonsocial_lvl2[i] <- read_meants("AIns_cat2_ugr_meants", s)
  out$nonsocial_lvl3[i] <- read_meants("AIns_cat3_ugr_meants", s)
  out$nonsocial_lvl4[i] <- read_meants("AIns_cat4_ugr_meants", s)
  out$nonsocial_lvl1_var[i] <- read_meants("AIns_var_cat1_ugr_meants", s)
  out$nonsocial_lvl2_var[i] <- read_meants("AIns_var_cat2_ugr_meants", s)
  out$nonsocial_lvl3_var[i] <- read_meants("AIns_var_cat3_ugr_meants", s)
  out$nonsocial_lvl4_var[i] <- read_meants("AIns_var_cat4_ugr_meants", s)
  out$social_lvl1[i] <- read_meants("AIns_cat5_ugr_meants", s)
  out$social_lvl2[i] <- read_meants("AIns_cat6_ugr_meants", s)
  out$social_lvl3[i] <- read_meants("AIns_cat7_ugr_meants", s)
  out$social_lvl4[i] <- read_meants("AIns_cat8_ugr_meants", s)
  out$social_lvl1_var[i] <- read_meants("AIns_var_cat5_ugr_meants", s)
  out$social_lvl2_var[i] <- read_meants("AIns_var_cat6_ugr_meants", s)
  out$social_lvl3_var[i] <- read_meants("AIns_var_cat7_ugr_meants", s)
  out$social_lvl4_var[i] <- read_meants("AIns_var_cat8_ugr_meants", s)
}

# inspect result
out





# --- after your loop finishes and `out` exists ---

library(tidyverse)

# helper to make mean +/- SE bar plot for a set of columns
plot_bar_means <- function(df, cols, title) {
  df %>%
    select(all_of(cols)) %>%
    pivot_longer(cols = everything(), names_to = "level", values_to = "value") %>%
    group_by(level) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      .groups = "drop"
    ) %>%
    mutate(level = factor(level, levels = cols)) %>%
    ggplot(aes(x = level, y = mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    labs(x = NULL, y = "Mean AIns meants (±SE)", title = title) +
    theme_classic(base_size = 14)
}

# column sets
nonsocial_cols <- c("nonsocial_lvl1","nonsocial_lvl2","nonsocial_lvl3","nonsocial_lvl4")
social_cols    <- c("social_lvl1","social_lvl2","social_lvl3","social_lvl4")

# plots
p_nonsocial <- plot_bar_means(out, nonsocial_cols, "Nonsocial (Levels 1–4)")
p_social    <- plot_bar_means(out, social_cols,    "Social (Levels 1–4)")

p_nonsocial
p_social




library(tidyverse)

make_invvar_weights <- function(var_vec, var_floor = 1e-6) {
  1 / pmax(var_vec, var_floor)
}

summarize_weighted_like_example <- function(df, mean_cols, var_cols, var_floor = 1e-6) {
  stopifnot(length(mean_cols) == length(var_cols))
  
  purrr::map2_dfr(mean_cols, var_cols, \(m, v) {
    y <- df[[m]]
    vv <- df[[v]]
    w <- make_invvar_weights(vv, var_floor = var_floor)
    
    ok <- is.finite(y) & is.finite(w) & (w > 0)
    y <- y[ok]; w <- w[ok]; vv <- vv[ok]
    
    tibble(
      level = m,
      n = length(y),
      
      # EXACTLY their weighted mean
      w_mean = if (length(y) > 0) sum(w * y) / sum(w) else NA_real_,
      
      # EXACTLY their SEM
      w_se = if (length(y) > 0) sqrt(1 / sum(w)) else NA_real_,
      
      # diagnostics to reveal the problem
      var_min = min(vv, na.rm = TRUE),
      var_p01 = as.numeric(quantile(vv, 0.01, na.rm = TRUE)),
      var_med = median(vv, na.rm = TRUE),
      var_p99 = as.numeric(quantile(vv, 0.99, na.rm = TRUE)),
      var_max = max(vv, na.rm = TRUE),
      
      w_sum = sum(w),
      w_max = max(w),
      w_p99 = as.numeric(quantile(w, 0.99, na.rm = TRUE))
    )
  }) %>%
    mutate(level = paste0("lvl", stringr::str_extract(level, "\\d+")))
}

plot_weighted_bars_like_example <- function(sum_df, title) {
  ggplot(sum_df, aes(x = level, y = w_mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = w_mean - w_se, ymax = w_mean + w_se), width = 0.2) +
    labs(x = NULL, y = "Weighted mean (±SEM), SEM = sqrt(1/sum(w))", title = title) +
    theme_classic(base_size = 14)
}

# columns
nonsocial_mean_cols <- c("nonsocial_lvl1","nonsocial_lvl2","nonsocial_lvl3","nonsocial_lvl4")
nonsocial_var_cols  <- c("nonsocial_lvl1_var","nonsocial_lvl2_var","nonsocial_lvl3_var","nonsocial_lvl4_var")

social_mean_cols <- c("social_lvl1","social_lvl2","social_lvl3","social_lvl4")
social_var_cols  <- c("social_lvl1_var","social_lvl2_var","social_lvl3_var","social_lvl4_var")

# summaries (EXACT example method)
sum_nonsocial_ex <- summarize_weighted_like_example(out, nonsocial_mean_cols, nonsocial_var_cols)
sum_social_ex    <- summarize_weighted_like_example(out, social_mean_cols, social_var_cols)

# print diagnostics
sum_nonsocial_ex
sum_social_ex

# plots
plot_weighted_bars_like_example(sum_nonsocial_ex, "Nonsocial (Levels 1–4): weighted (example method)")
plot_weighted_bars_like_example(sum_social_ex,    "Social (Levels 1–4): weighted (example method)")
