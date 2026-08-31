#!/usr/bin/env Rscript

# Descriptive voxelwise Trust Game--Ultimatum Game correspondence in an
# independent anatomical left-insula ROI. This script only reads the final
# subject-level images routed into the existing N=132 L3 analyses. It does not
# rerun or modify any FEAT model and never extracts from an L3 statistical map.

options(stringsAsFactors = FALSE, warn = 1)

TG_EXPECTED_COPE <- 10L
UG_EXPECTED_COPE <- 11L
GRID_TOLERANCE <- 1e-5

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Run this analysis with Rscript.", call. = FALSE)
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

usage <- function(default_root) {
  cat(paste0(
    "Usage: Rscript code/voxelwise_cross_task_insula.R [options]\n\n",
    "Options:\n",
    "  --project-root PATH       rf1-betrayal root [", default_root, "]\n",
    "  --output-dir PATH         statistical output directory [code/voxelwise_cross_task_insula_output]\n",
    "  --fsldir PATH             FSL root [$FSLDIR]\n",
    "  --atlas-file PATH         explicit cortical maxprob-thr25 atlas (optional)\n",
    "  --atlas-xml PATH          explicit Harvard-Oxford cortical XML (optional)\n",
    "  --help                    show this message\n"
  ))
}

parse_args <- function() {
  default_root <- normalizePath(file.path(dirname(script_path()), ".."), mustWork = TRUE)
  result <- list(
    project_root = default_root,
    output_dir = NULL,
    fsldir = Sys.getenv("FSLDIR", unset = NA_character_),
    atlas_file = NULL,
    atlas_xml = NULL
  )
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      usage(default_root)
      quit(status = 0L)
    }
    key <- switch(
      arg,
      "--project-root" = "project_root",
      "--output-dir" = "output_dir",
      "--fsldir" = "fsldir",
      "--atlas-file" = "atlas_file",
      "--atlas-xml" = "atlas_xml",
      NULL
    )
    if (is.null(key)) stopf("Unknown argument: %s", arg)
    if (i == length(args)) stopf("Missing value after %s", arg)
    result[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  result$project_root <- normalizePath(result$project_root, mustWork = TRUE)
  if (is.null(result$output_dir)) {
    result$output_dir <- file.path(
      result$project_root, "code", "voxelwise_cross_task_insula_output"
    )
  }
  result
}

require_packages <- function() {
  packages <- c("RNifti", "xml2", "ggplot2", "jsonlite")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stopf(
      "Missing R package(s): %s. Install them before running this analysis.",
      paste(missing, collapse = ", ")
    )
  }
  invisible(packages)
}

read_subjects <- function(path) {
  if (!file.exists(path)) stopf("Subject list not found: %s", path)
  subjects <- trimws(readLines(path, warn = FALSE))
  subjects <- sub("^sub-", "", subjects)
  subjects <- subjects[nzchar(subjects)]
  if (!length(subjects)) stopf("Subject list is empty: %s", path)
  if (anyDuplicated(subjects)) stopf("Subject list contains duplicate IDs: %s", path)
  subjects
}

regex_rows <- function(lines, pattern) {
  hits <- regexec(pattern, lines, perl = TRUE)
  values <- regmatches(lines, hits)
  values[lengths(values) > 0L]
}

parse_contrast_names <- function(path) {
  rows <- regex_rows(
    readLines(path, warn = FALSE),
    '^set fmri\\(conname_real\\.([0-9]+)\\) "([^"]+)"$'
  )
  setNames(vapply(rows, `[[`, character(1), 3L), vapply(rows, `[[`, character(1), 2L))
}

parse_ev_titles <- function(path) {
  rows <- regex_rows(
    readLines(path, warn = FALSE),
    '^set fmri\\(evtitle([0-9]+)\\) "([^"]+)"$'
  )
  setNames(vapply(rows, `[[`, character(1), 3L), vapply(rows, `[[`, character(1), 2L))
}

parse_weights <- function(path, cope) {
  rows <- regex_rows(
    readLines(path, warn = FALSE),
    sprintf('^set fmri\\(con_real%d\\.([0-9]+)\\)[[:space:]]+([-+0-9.eE]+)$', cope)
  )
  setNames(as.numeric(vapply(rows, `[[`, character(1), 3L)), vapply(rows, `[[`, character(1), 2L))
}

find_named_contrast <- function(contrast_names, exact_name, path) {
  found <- as.integer(names(contrast_names)[unname(contrast_names) == exact_name])
  if (length(found) != 1L) {
    stopf("Expected one contrast named '%s' in %s; found: %s", exact_name, path, paste(found, collapse = ", "))
  }
  found
}

validate_contrasts <- function(project_root) {
  tg_fsf <- file.path(project_root, "templates", "L1_task-trust_model-01_type-act.fsf")
  ug_fsf <- file.path(project_root, "templates", "L1_task-ugr_model-3_type-act.fsf")
  ug_generator <- file.path(project_root, "code", "a4_model-3.py")
  required <- c(tg_fsf, ug_fsf, ug_generator)
  if (any(!file.exists(required))) stopf("Missing contrast-definition file: %s", required[!file.exists(required)][1])

  tg_cope <- find_named_contrast(parse_contrast_names(tg_fsf), "rec-def", tg_fsf)
  ug_cope <- find_named_contrast(
    parse_contrast_names(ug_fsf), "offer (un)fairness (pmod)", ug_fsf
  )
  if (tg_cope != TG_EXPECTED_COPE || ug_cope != UG_EXPECTED_COPE) {
    stopf(
      "Contrast numbering changed: expected TG cope 10 and UGR cope 11; found %d and %d.",
      tg_cope, ug_cope
    )
  }

  tg_titles <- parse_ev_titles(tg_fsf)
  tg_weights <- parse_weights(tg_fsf, tg_cope)
  defect_evs <- as.integer(names(tg_titles)[grepl("_def$", tg_titles)])
  recip_evs <- as.integer(names(tg_titles)[grepl("_rec$", tg_titles)])
  if (!identical(sort(defect_evs), c(4L, 6L, 8L)) ||
      !identical(sort(recip_evs), c(5L, 7L, 9L)) ||
      !all(abs(tg_weights[as.character(defect_evs)] + 1) < 1e-12) ||
      !all(abs(tg_weights[as.character(recip_evs)] - 1) < 1e-12)) {
    stop("TG cope 10 is no longer reciprocated > unreciprocated; review sign handling.", call. = FALSE)
  }

  ug_titles <- parse_ev_titles(ug_fsf)
  ug_weights <- parse_weights(ug_fsf, ug_cope)
  offer_pmods <- as.integer(names(ug_titles)[grepl("_pmod$", ug_titles) & !grepl("^rt_", ug_titles)])
  if (!identical(offer_pmods, c(2L, 4L, 6L, 8L)) ||
      !all(abs(ug_weights[as.character(offer_pmods)] - 1) < 1e-12)) {
    stop("UGR cope 11 no longer positively weights all four offer pmods.", call. = FALSE)
  }
  generator <- paste(readLines(ug_generator, warn = FALSE), collapse = "\n")
  checks <- c(
    'df["offer_amount"] = df["L_Option"] + df["R_Option"]',
    'group["demeaned_offer"] = group["offer_amount"] - group["offer_amount"].mean()',
    '"demeaned_offer"'
  )
  if (!all(vapply(checks, grepl, logical(1), x = generator, fixed = TRUE))) {
    stopf("Could not verify positive demeaned-offer pmods in %s", ug_generator)
  }

  list(
    tg_cope = tg_cope,
    tg_name = "rec-def",
    tg_stored_direction = "reciprocated > unreciprocated",
    tg_applied_multiplier = 1,
    tg_l1_fsf = tg_fsf,
    ug_cope = ug_cope,
    ug_name = "offer (un)fairness (pmod)",
    ug_stored_direction = "fairness (positive demeaned offer amount)",
    ug_applied_multiplier = 1,
    ug_l1_fsf = ug_fsf,
    ug_pmod_generator = ug_generator
  )
}

remap_template_path <- function(path, project_root) {
  if (!grepl("/derivatives/", path, fixed = TRUE)) {
    stopf("L3 input is not under a derivatives directory: %s", path)
  }
  relative <- sub("^.*?/derivatives/", "", path)
  file.path(project_root, "derivatives", relative)
}

parse_l3_routes <- function(template, cope, project_root) {
  if (!file.exists(template)) stopf("L3 template not found: %s", template)
  rows <- regex_rows(
    readLines(template, warn = FALSE),
    '^set feat_files\\(([0-9]+)\\) "([^"]+)"$'
  )
  if (!length(rows)) stopf("No feat_files entries found in %s", template)
  indices <- as.integer(vapply(rows, `[[`, character(1), 2L))
  if (!identical(indices, seq_along(indices))) stopf("Non-contiguous feat_files entries in %s", template)
  raw_paths <- gsub("COPENUM", as.character(cope), vapply(rows, `[[`, character(1), 3L), fixed = TRUE)
  subjects <- sub("^.*/sub-([^/]+)/.*$", "\\1", raw_paths)
  if (any(subjects == raw_paths)) stopf("Could not parse subject IDs from %s", template)
  cope_paths <- vapply(raw_paths, remap_template_path, character(1), project_root = project_root)
  zstat_paths <- sub("/cope([0-9]+)\\.nii\\.gz$", "/zstat\\1.nii.gz", cope_paths)
  if (any(zstat_paths == cope_paths)) stopf("Unexpected COPE filename in %s", template)
  mask_paths <- file.path(dirname(dirname(cope_paths)), "mask.nii.gz")
  data.frame(
    subject = subjects,
    template_path = raw_paths,
    cope = unname(cope_paths),
    zstat = unname(zstat_paths),
    mask = unname(mask_paths),
    level = ifelse(grepl("/L1_", raw_paths, fixed = TRUE), "L1", "L2"),
    run = ifelse(
      grepl("_run-[0-9]+_", raw_paths),
      sub("^.*_run-([0-9]+)_.*$", "\\1", raw_paths),
      "both"
    )
  )
}

assert_route_agreement <- function(subjects, primary, alternate, task) {
  if (!identical(primary$subject, subjects)) {
    stopf("%s L3 template order does not exactly match sublist_n132.txt.", task)
  }
  if (!identical(primary[, c("subject", "template_path")], alternate[, c("subject", "template_path")])) {
    stopf("%s full and ones N=132 L3 templates use different inputs.", task)
  }
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

require_complete_sample <- function(subjects, tg, ug, output_dir) {
  required_labels <- c("TG COPE", "TG ZSTAT", "TG FEAT mask", "UGR COPE", "UGR ZSTAT", "UGR FEAT mask")
  included <- logical(length(subjects))
  excluded_rows <- vector("list", length(subjects))
  for (i in seq_along(subjects)) {
    paths <- c(tg$cope[i], tg$zstat[i], tg$mask[i], ug$cope[i], ug$zstat[i], ug$mask[i])
    missing <- !file.exists(paths)
    if (any(missing)) {
      reason <- paste(sprintf("%s missing: %s", required_labels[missing], paths[missing]), collapse = "; ")
      excluded_rows[[i]] <- data.frame(subject = subjects[i], reason = reason)
      message(sprintf("EXCLUDED sub-%s: %s", subjects[i], reason))
    } else {
      included[i] <- TRUE
    }
  }
  excluded <- if (any(!included)) do.call(rbind, excluded_rows[!included]) else data.frame(subject = character(), reason = character())
  included_audit <- data.frame(
    subject = subjects[included],
    tg_level = tg$level[included], tg_run = tg$run[included],
    tg_cope = tg$cope[included], tg_zstat = tg$zstat[included], tg_mask = tg$mask[included],
    ug_level = ug$level[included], ug_run = ug$run[included],
    ug_cope = ug$cope[included], ug_zstat = ug$zstat[included], ug_mask = ug$mask[included]
  )
  write_tsv(excluded, file.path(output_dir, "excluded_subjects.tsv"))
  write_tsv(included_audit, file.path(output_dir, "included_subjects.tsv"))
  if (nrow(excluded)) {
    stopf(
      paste0(
        "%d of the required 132 subjects lack one or more inputs. No analysis was run and ",
        "the sample was not reduced; see %s"
      ),
      nrow(excluded), file.path(output_dir, "excluded_subjects.tsv")
    )
  }
  if (!all(included) || length(subjects) != 132L) {
    stop("Internal sample validation failed: the analysis requires exactly all 132 subjects.", call. = FALSE)
  }
  list(tg = tg, ug = ug, excluded = excluded)
}

load_image <- function(path, internal = FALSE) {
  image <- RNifti::readNifti(path, internal = internal)
  if (length(dim(image)) != 3L) stopf("Expected 3D image; found %s: %s", paste(dim(image), collapse = "x"), path)
  image
}

affine_of <- function(image) {
  matrix(as.numeric(RNifti::xform(image, useQuaternionFirst = FALSE)), nrow = 4L, ncol = 4L)
}

same_grid <- function(a, b) {
  identical(as.integer(dim(a)), as.integer(dim(b))) &&
    max(abs(affine_of(a) - affine_of(b))) <= GRID_TOLERANCE
}

validate_subject_grids <- function(tg, ug) {
  reference_path <- tg$cope[1]
  reference <- load_image(reference_path)
  for (i in seq_len(nrow(tg))) {
    paths <- c(tg$cope[i], tg$zstat[i], tg$mask[i], ug$cope[i], ug$zstat[i], ug$mask[i])
    labels <- c("TG COPE", "TG ZSTAT", "TG mask", "UGR COPE", "UGR ZSTAT", "UGR mask")
    for (j in seq_along(paths)) {
      image <- load_image(paths[j], internal = TRUE)
      if (!same_grid(reference, image)) {
        stopf(
          paste0(
            "Grid mismatch for sub-%s %s: %s\nReference %s shape=%s affine=\n%s\n",
            "Image shape=%s affine=\n%s\nSubject statistical images are never resampled."
          ),
          tg$subject[i], labels[j], paths[j], reference_path,
          paste(dim(reference), collapse = "x"), paste(capture.output(print(affine_of(reference))), collapse = "\n"),
          paste(dim(image), collapse = "x"), paste(capture.output(print(affine_of(image))), collapse = "\n")
        )
      }
    }
  }
  reference
}

discover_xml <- function(fsldir, explicit = NULL) {
  if (!is.null(explicit)) {
    if (!file.exists(explicit)) stopf("Atlas XML not found: %s", explicit)
    return(normalizePath(explicit))
  }
  atlas_dir <- file.path(fsldir, "data", "atlases")
  preferred <- file.path(atlas_dir, "HarvardOxford-Cortical.xml")
  if (file.exists(preferred)) return(normalizePath(preferred))
  candidates <- list.files(atlas_dir, pattern = "HarvardOxford.*Cortical.*\\.xml$", recursive = TRUE, full.names = TRUE)
  if (length(candidates) != 1L) {
    stopf("Could not uniquely locate Harvard-Oxford cortical XML under %s; found %d", atlas_dir, length(candidates))
  }
  normalizePath(candidates)
}

atlas_candidates <- function(fsldir, explicit = NULL) {
  if (!is.null(explicit)) {
    if (!file.exists(explicit)) stopf("Atlas NIfTI not found: %s", explicit)
    name <- basename(explicit)
    if (!grepl("HarvardOxford-cort-maxprob-thr25-", name, fixed = TRUE)) {
      stopf("--atlas-file is not a cortical maxprob-thr25 atlas: %s", name)
    }
    return(normalizePath(explicit))
  }
  atlas_dir <- file.path(fsldir, "data", "atlases")
  candidates <- list.files(
    atlas_dir,
    pattern = "^HarvardOxford-cort-maxprob-thr25-.*\\.nii(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(candidates)) stopf("No cortical Harvard-Oxford maxprob-thr25 atlas found under %s", atlas_dir)
  normalizePath(sort(candidates))
}

choose_atlas <- function(candidates, reference) {
  # Inspect candidate headers without retaining every atlas array in memory.
  images <- lapply(candidates, load_image, internal = TRUE)
  exact <- which(vapply(images, same_grid, logical(1), b = reference))
  if (length(exact)) {
    chosen <- exact[1]
    rule <- "exact statistical-grid match"
  } else {
    volumes <- vapply(images, function(x) abs(det(affine_of(x)[1:3, 1:3, drop = FALSE])), numeric(1))
    chosen <- order(volumes, candidates)[1]
    rule <- "no exact grid match; discovered atlas with smallest voxel volume"
  }
  list(path = candidates[chosen], image = load_image(candidates[chosen]), rule = rule)
}

insula_label_from_xml <- function(path) {
  doc <- xml2::read_xml(path)
  labels <- xml2::xml_find_all(doc, ".//label[normalize-space(text())='Insular Cortex']")
  if (length(labels) != 1L) stopf("Expected one XML label named 'Insular Cortex' in %s; found %d", path, length(labels))
  index <- as.integer(xml2::xml_attr(labels, "index"))
  if (is.na(index)) stopf("Insular Cortex label has no integer index in %s", path)
  # FSL XML indices are zero-based; max-probability atlas value 0 is background.
  list(name = trimws(xml2::xml_text(labels)), xml_index = index, nifti_value = index + 1L)
}

voxel_to_world <- function(indices_one_based, affine) {
  homogeneous <- cbind(indices_one_based - 1, 1)
  transformed <- homogeneous %*% t(affine)
  transformed[, 1:3, drop = FALSE]
}

nearest_resample_mask <- function(source_mask, source_affine, target_dim, target_affine, chunk_size = 250000L) {
  output <- array(FALSE, dim = target_dim)
  inverse_source <- solve(source_affine)
  total <- prod(target_dim)
  starts <- seq.int(1L, total, by = chunk_size)
  for (start in starts) {
    linear <- start:min(start + chunk_size - 1L, total)
    target_index <- arrayInd(linear, .dim = target_dim)
    world <- voxel_to_world(target_index, target_affine)
    source_zero <- cbind(world, 1) %*% t(inverse_source)
    source_index <- round(source_zero[, 1:3, drop = FALSE]) + 1L
    valid <- apply(source_index >= 1L & sweep(source_index, 2L, dim(source_mask), `<=`), 1L, all)
    values <- logical(length(linear))
    if (any(valid)) values[valid] <- source_mask[source_index[valid, , drop = FALSE]]
    output[linear] <- values
  }
  output
}

make_left_insula_mask <- function(atlas, label_value, reference) {
  atlas_data <- as.array(atlas)
  assigned_index <- which(abs(atlas_data - label_value) < 1e-6, arr.ind = TRUE)
  if (!nrow(assigned_index)) stopf("Atlas label %d is empty.", label_value)
  xyz <- voxel_to_world(assigned_index, affine_of(atlas))
  left_index <- assigned_index[xyz[, 1] < 0, , drop = FALSE]
  native <- array(FALSE, dim = dim(atlas_data))
  native[left_index] <- TRUE
  native_count <- sum(native)
  if (!native_count) stopf("Atlas label %d produced an empty left-hemisphere mask.", label_value)
  resampled <- !same_grid(atlas, reference)
  reference_mask <- if (resampled) {
    nearest_resample_mask(native, affine_of(atlas), dim(reference), affine_of(reference))
  } else native
  if (!sum(reference_mask)) stop("Left-insula mask is empty on the statistical reference grid.", call. = FALSE)
  list(mask = reference_mask, native_count = native_count, reference_count = sum(reference_mask), resampled = resampled)
}

coverage_intersection <- function(tg, ug, reference_dim) {
  tg_coverage <- array(TRUE, dim = reference_dim)
  ug_coverage <- array(TRUE, dim = reference_dim)
  for (i in seq_len(nrow(tg))) {
    tg_mask <- as.array(load_image(tg$mask[i]))
    ug_mask <- as.array(load_image(ug$mask[i]))
    tg_coverage <- tg_coverage & is.finite(tg_mask) & tg_mask > 0
    ug_coverage <- ug_coverage & is.finite(ug_mask) & ug_mask > 0
  }
  list(tg = tg_coverage, ug = ug_coverage)
}

write_nifti <- function(data, path, reference, datatype) {
  RNifti::writeNifti(data, path, template = reference, datatype = datatype)
  invisible(path)
}

running_stats <- function(size) list(n = 0L, mean = numeric(size), m2 = numeric(size))

update_stats <- function(stat, values) {
  stat$n <- stat$n + 1L
  delta <- values - stat$mean
  stat$mean <- stat$mean + delta / stat$n
  stat$m2 <- stat$m2 + delta * (values - stat$mean)
  stat
}

finish_stats <- function(stat) {
  if (stat$n < 2L) stop("At least two included subjects are required.", call. = FALSE)
  sd <- sqrt(pmax(stat$m2 / (stat$n - 1L), 0))
  list(mean = stat$mean, sd = sd, sem = sd / sqrt(stat$n))
}

extract_summaries <- function(tg, ug, final_mask) {
  size <- sum(final_mask)
  accum <- setNames(replicate(4L, running_stats(size), simplify = FALSE), c("cope_tg", "cope_ug", "zstat_tg", "zstat_ug"))
  for (i in seq_len(nrow(tg))) {
    values <- list(
      cope_tg = as.numeric(load_image(tg$cope[i])[final_mask]),
      cope_ug = as.numeric(load_image(ug$cope[i])[final_mask]),
      zstat_tg = as.numeric(load_image(tg$zstat[i])[final_mask]),
      zstat_ug = as.numeric(load_image(ug$zstat[i])[final_mask])
    )
    for (name in names(values)) {
      if (any(!is.finite(values[[name]]))) {
        stopf("sub-%s has %d non-finite %s values inside the final ROI.", tg$subject[i], sum(!is.finite(values[[name]])), name)
      }
      accum[[name]] <- update_stats(accum[[name]], values[[name]])
    }
  }
  list(
    cope = list(tg = finish_stats(accum$cope_tg), ug = finish_stats(accum$cope_ug)),
    zstat = list(tg = finish_stats(accum$zstat_tg), ug = finish_stats(accum$zstat_ug))
  )
}

fit_odr <- function(x, y) {
  if (sd(x) < .Machine$double.eps || sd(y) < .Machine$double.eps) {
    stop("Cannot fit ODR because a voxelwise mean has zero variance.", call. = FALSE)
  }
  centered <- scale(cbind(x, y), center = TRUE, scale = FALSE)
  eig <- eigen(crossprod(centered), symmetric = TRUE)
  direction <- eig$vectors[, 1L]
  if (direction[1] < 0) direction <- -direction
  if (abs(direction[1]) < 1e-12) stop("ODR solution is vertical and cannot be expressed as y = a + bx.", call. = FALSE)
  slope <- direction[2] / direction[1]
  intercept <- mean(y) - slope * mean(x)
  residual <- (y - (intercept + slope * x)) / sqrt(1 + slope^2)
  list(
    pearson_r = unname(cor(x, y, method = "pearson")),
    odr_intercept = unname(intercept),
    odr_slope = unname(slope),
    signed_orthogonal_residual = unname(residual),
    method = "unweighted orthogonal distance regression (total least squares)"
  )
}

write_voxel_table <- function(path, quantity, indices, coordinates, summary, fit) {
  prefix <- if (quantity == "cope") "cope" else "subject_zstat"
  result <- data.frame(
    i = indices[, 1] - 1L,
    j = indices[, 2] - 1L,
    k = indices[, 3] - 1L,
    mni_x_mm = coordinates[, 1],
    mni_y_mm = coordinates[, 2],
    mni_z_mm = coordinates[, 3]
  )
  result[[sprintf("tg_mean_%s_recip_gt_nonrecip", prefix)]] <- summary$tg$mean
  result[[sprintf("tg_sd_%s_recip_gt_nonrecip", prefix)]] <- summary$tg$sd
  result[[sprintf("tg_sem_%s_recip_gt_nonrecip", prefix)]] <- summary$tg$sem
  result[[sprintf("ug_mean_%s_fairness_pmod", prefix)]] <- summary$ug$mean
  result[[sprintf("ug_sd_%s_fairness_pmod", prefix)]] <- summary$ug$sd
  result[[sprintf("ug_sem_%s_fairness_pmod", prefix)]] <- summary$ug$sem
  result$signed_odr_residual <- fit$signed_orthogonal_residual
  write_tsv(result, path)
}

make_plot <- function(stem, quantity, x, y, fit, n_subjects) {
  label <- if (quantity == "cope") "mean COPE" else "mean subject-level Z-statistic"
  data <- data.frame(tg = x, ug = y)
  annotation <- sprintf(
    "N = %d subjects; V = %d voxels\nPearson r = %.3f\nODR: y = %.3g + %.3gx",
    n_subjects, length(x), fit$pearson_r, fit$odr_intercept, fit$odr_slope
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = tg, y = ug)) +
    ggplot2::geom_point(shape = 21, size = 1.8, alpha = 0.58, fill = "#276FBF", colour = "white", stroke = 0.2) +
    ggplot2::geom_abline(intercept = fit$odr_intercept, slope = fit$odr_slope, colour = "#B23A48", linewidth = 0.9) +
    ggplot2::annotate("label", x = -Inf, y = Inf, label = annotation, hjust = -0.05, vjust = 1.05, size = 3.2, label.size = 0.2) +
    ggplot2::labs(
      title = "Voxelwise cross-task correspondence in anatomical left insula",
      x = sprintf("Trust Game: reciprocated > nonreciprocated\n(%s)", label),
      y = sprintf("Ultimatum Game: fairness parametric modulation\n(%s)", label),
      caption = "Descriptive voxelwise correspondence; spatial dependence precludes treating voxels as independent."
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.caption = ggplot2::element_text(size = 8, colour = "grey35"),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3)
    )
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = 7.1, height = 6.0, units = "in")
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = 7.1, height = 6.0, units = "in", dpi = 300)
}

main <- function() {
  args <- parse_args()
  require_packages()
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(args$output_dir, mustWork = TRUE)
  masks_dir <- file.path(args$project_root, "masks")
  if (!dir.exists(masks_dir)) stopf("Repository masks directory not found: %s", masks_dir)
  masks_dir <- normalizePath(masks_dir, mustWork = TRUE)
  if (is.na(args$fsldir) || !nzchar(args$fsldir)) stop("FSLDIR is not set. Export FSLDIR or pass --fsldir.", call. = FALSE)
  fsldir <- normalizePath(args$fsldir, mustWork = TRUE)

  message("Validating contrast definitions and L3 subject-level routes...")
  contrast <- validate_contrasts(args$project_root)
  subjects_path <- file.path(args$project_root, "code", "sublist_n132.txt")
  subjects <- read_subjects(subjects_path)
  if (length(subjects) != 132L) stopf("Expected 132 subjects; found %d in %s", length(subjects), subjects_path)

  templates <- file.path(args$project_root, "templates")
  template_paths <- list(
    tg_full = file.path(templates, "L3_task-trust_model-01_type-act_group-full_n132_flame1.fsf"),
    tg_ones = file.path(templates, "L3_task-trust_model-01_type-act_group-ones_n132_flame1.fsf"),
    ug_full = file.path(templates, "L3_task-ugr_model-3_type-act_group-full_n132_flame1.fsf"),
    ug_ones = file.path(templates, "L3_task-ugr_model-3_type-act_group-ones_n132_flame1.fsf")
  )
  tg <- parse_l3_routes(template_paths$tg_full, contrast$tg_cope, args$project_root)
  ug <- parse_l3_routes(template_paths$ug_full, contrast$ug_cope, args$project_root)
  tg_ones <- parse_l3_routes(template_paths$tg_ones, contrast$tg_cope, args$project_root)
  ug_ones <- parse_l3_routes(template_paths$ug_ones, contrast$ug_cope, args$project_root)
  assert_route_agreement(subjects, tg, tg_ones, "TG")
  assert_route_agreement(subjects, ug, ug_ones, "UGR")

  cases <- require_complete_sample(subjects, tg, ug, output_dir)
  tg <- cases$tg
  ug <- cases$ug
  message(sprintf("Final subject N: %d (all subjects in sublist_n132.txt are required)", nrow(tg)))
  reference <- validate_subject_grids(tg, ug)
  reference_path <- tg$cope[1]
  message(sprintf("Statistical reference: %s", reference_path))
  message(sprintf("Reference shape: %s", paste(dim(reference), collapse = " x ")))
  message("Reference affine:\n", paste(capture.output(print(affine_of(reference))), collapse = "\n"))

  xml_path <- discover_xml(fsldir, args$atlas_xml)
  atlas_choice <- choose_atlas(atlas_candidates(fsldir, args$atlas_file), reference)
  label <- insula_label_from_xml(xml_path)
  insula <- make_left_insula_mask(atlas_choice$image, label$nifti_value, reference)
  message(sprintf("Harvard-Oxford atlas: %s", atlas_choice$path))
  message(sprintf("Atlas selection: %s", atlas_choice$rule))
  message(sprintf("Atlas XML: %s", xml_path))
  message(sprintf("Atlas XML label: '%s', XML index=%d, NIfTI integer=%d", label$name, label$xml_index, label$nifti_value))
  message(sprintf("Initial whole-left-insula voxels (atlas grid): %d", insula$native_count))
  if (insula$resampled) message(sprintf("Whole-left-insula voxels after nearest-neighbor mask resampling: %d", insula$reference_count))

  coverage <- coverage_intersection(tg, ug, dim(reference))
  final_mask <- insula$mask & coverage$tg & coverage$ug
  final_voxels <- sum(final_mask)
  if (!final_voxels) stop("Anatomical left-insula x TG coverage x UGR coverage mask is empty.", call. = FALSE)
  tg_covered <- sum(insula$mask & coverage$tg)
  ug_covered <- sum(insula$mask & coverage$ug)
  message(sprintf("Left-insula voxels with all-subject TG coverage: %d", tg_covered))
  message(sprintf("Left-insula voxels with all-subject UGR coverage: %d", ug_covered))
  message(sprintf("Final dual-task coverage voxel count: %d", final_voxels))

  anatomical_path <- file.path(masks_dir, "left_insula_maxprob-thr25_anatomical_refgrid.nii.gz")
  final_mask_path <- file.path(masks_dir, "left_insula_maxprob-thr25_dualtask_coverage_mask.nii.gz")
  write_nifti(insula$mask * 1L, anatomical_path, reference, "uint8")
  write_nifti(final_mask * 1L, final_mask_path, reference, "uint8")

  indices <- which(final_mask, arr.ind = TRUE)
  coordinates <- voxel_to_world(indices, affine_of(reference))
  summaries <- extract_summaries(tg, ug, final_mask)
  fits <- list()
  output_files <- list(
    anatomical_mask_reference_grid = anatomical_path,
    final_dualtask_coverage_mask = final_mask_path,
    included_subjects = file.path(output_dir, "included_subjects.tsv"),
    excluded_subjects = file.path(output_dir, "excluded_subjects.tsv")
  )
  for (quantity in c("cope", "zstat")) {
    summary <- summaries[[quantity]]
    fit <- fit_odr(summary$tg$mean, summary$ug$mean)
    fits[[quantity]] <- fit
    table_path <- file.path(output_dir, sprintf("voxelwise_cross_task_insula_%s.tsv", quantity))
    plot_stem <- file.path(output_dir, sprintf("voxelwise_cross_task_insula_%s_scatter", quantity))
    residual_path <- file.path(output_dir, sprintf("voxelwise_cross_task_insula_%s_signed_odr_residual.nii.gz", quantity))
    write_voxel_table(table_path, quantity, indices, coordinates, summary, fit)
    make_plot(plot_stem, quantity, summary$tg$mean, summary$ug$mean, fit, nrow(tg))
    residual_map <- array(0, dim = dim(reference))
    residual_map[final_mask] <- fit$signed_orthogonal_residual
    write_nifti(residual_map, residual_path, reference, "float")
    output_files[[paste0(quantity, "_voxel_table")]] <- table_path
    output_files[[paste0(quantity, "_scatter_pdf")]] <- paste0(plot_stem, ".pdf")
    output_files[[paste0(quantity, "_scatter_png")]] <- paste0(plot_stem, ".png")
    output_files[[paste0(quantity, "_signed_odr_residual_map")]] <- residual_path
    display <- if (quantity == "cope") "COPE" else "mean subject-level ZSTAT"
    message(sprintf(
      "%s: Pearson r=%.6f; ODR slope=%.6f; intercept=%.6f",
      display, fit$pearson_r, fit$odr_slope, fit$odr_intercept
    ))
  }

  fit_metadata <- lapply(fits, function(x) x[setdiff(names(x), "signed_orthogonal_residual")])
  metadata_path <- file.path(output_dir, "analysis_metadata.json")
  output_files$metadata <- metadata_path
  metadata <- list(
    analysis = "descriptive voxelwise cross-task correspondence in anatomical left insula",
    implementation = "R",
    package_versions = as.list(vapply(c("RNifti", "xml2", "ggplot2", "jsonlite"), function(x) as.character(utils::packageVersion(x)), character(1))),
    project_root = args$project_root,
    subject_list = subjects_path,
    starting_n = length(subjects), included_n = nrow(tg), excluded_n = nrow(cases$excluded),
    contrast_validation = contrast,
    l3_templates = template_paths,
    reference_image = reference_path,
    reference_shape = as.integer(dim(reference)),
    reference_affine = unname(affine_of(reference)),
    fsl_dir = fsldir,
    harvard_oxford = list(
      atlas_file = atlas_choice$path,
      atlas_selection_rule = atlas_choice$rule,
      xml_file = xml_path,
      xml_label = label$name,
      xml_zero_based_index = label$xml_index,
      nifti_integer_value = label$nifti_value,
      left_hemisphere_rule = "voxel center MNI x < 0 mm",
      initial_left_insula_voxels_atlas_grid = insula$native_count,
      left_insula_voxels_reference_grid = insula$reference_count,
      mask_resampled_nearest_neighbor = insula$resampled
    ),
    coverage = list(
      criterion = "intersection of every included subject's final TG and UGR FEAT masks",
      left_insula_with_all_subject_tg_coverage = tg_covered,
      left_insula_with_all_subject_ugr_coverage = ug_covered,
      final_voxels = final_voxels
    ),
    voxel_index_convention = "i/j/k are zero-based NIfTI indices; MNI coordinates are millimetres",
    correspondence = fit_metadata,
    signed_residual_definition = paste(
      "(UG - (ODR intercept + ODR slope * TG)) / sqrt(1 + slope^2); positive means a",
      "stronger UG/fairness response than predicted from the TG reciprocity response, while",
      "negative means a weaker UG/fairness response than predicted. Residual maps are zero",
      "outside the saved final coverage mask."
    ),
    inference_note = paste(
      "Pearson r and ODR are descriptive. No conventional voxelwise Pearson p-value is computed",
      "because neighboring voxels are spatially dependent."
    ),
    outputs = output_files
  )
  jsonlite::write_json(metadata, metadata_path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  message(sprintf("Statistical outputs written to: %s", output_dir))
  message(sprintf("Masks written to: %s", masks_dir))
  invisible(metadata)
}

tryCatch(
  main(),
  error = function(error) {
    message("ERROR: ", conditionMessage(error))
    quit(status = 1L)
  }
)
