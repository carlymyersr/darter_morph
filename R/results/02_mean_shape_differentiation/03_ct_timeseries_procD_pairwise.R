# ============================================================
# Scripts/CT3_stats_models.R
# CT 3-TIMEPOINT STATISTICAL MODELS (RAW + RESIDUAL SHAPES)
#
# Canonical objects created upstream:
#   - R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R
#   - R/methods/01_specimen_sampling_study_design/01_build_metadata.R
#   - R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R
#
# Expected objects from R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R:
#   - gdf_ct3
#   - coords_ct3
#   - coords_resid_ct3
#
# Grouping variable used here:
#   - gdf_ct3$group (CT_1950, CT_1956, CT_1970)
# ============================================================

# ---------------------------
# Load canonical project objects
# ---------------------------
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R")


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(geomorph)
  library(RRPP)
  library(scales)
})


# ============================================================
# 0a) CT3 object aliasing (handles different naming conventions)
# ============================================================

# Common naming conventions I've seen in your earlier code:
#   gdf_3grp / coords_3grp / coords_resid_3grp
#   gdf_ct3  / coords_ct3  / coords_resid_ct3
#   gdf_CT3  / coords_CT3  / coords_resid_CT3

if (!exists("gdf_ct3")) {
  if (exists("gdf_3grp"))      gdf_ct3 <- gdf_3grp
  else if (exists("gdf_CT3"))  gdf_ct3 <- gdf_CT3
  else if (exists("gdf_ct_3grp")) gdf_ct3 <- gdf_ct_3grp
}

if (!exists("coords_ct3")) {
  if (exists("coords_3grp"))      coords_ct3 <- coords_3grp
  else if (exists("coords_CT3"))  coords_ct3 <- coords_CT3
  else if (exists("coords_ct_3grp")) coords_ct3 <- coords_ct_3grp
}

if (!exists("coords_resid_ct3")) {
  if (exists("coords_resid_3grp"))      coords_resid_ct3 <- coords_resid_3grp
  else if (exists("coords_resid_CT3"))  coords_resid_ct3 <- coords_resid_CT3
  else if (exists("coords_resid_ct_3grp")) coords_resid_ct3 <- coords_resid_ct_3grp
}

# If still missing, print the most likely candidates to help you pick the right names
if (!exists("gdf_ct3") || !exists("coords_ct3") || !exists("coords_resid_ct3")) {
  message("\nCT3 subset script ran, but expected object names not found.")
  message("Here are candidate objects in the environment:\n")
  
  objs <- ls(envir = .GlobalEnv)
  
  # likely metadata frames
  cand_df <- objs[vapply(objs, function(nm) {
    inherits(get(nm, envir = .GlobalEnv), "data.frame")
  }, logical(1))]
  
  # likely coords arrays
  cand_arr <- objs[vapply(objs, function(nm) {
    x <- get(nm, envir = .GlobalEnv)
    is.array(x) && length(dim(x)) == 3
  }, logical(1))]
  
  message("Data.frames (first 30): ", paste(head(cand_df, 30), collapse = ", "))
  message("3D arrays (first 30):   ", paste(head(cand_arr, 30), collapse = ", "))
  
  stop("\nPlease set gdf_ct3/coords_ct3/coords_resid_ct3 to the correct objects (see messages above).")
}





# ============================================================
# 0) Output directories + run ID + manifest
# ============================================================

run_id   <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_base <- file.path("Outputs", paste0("CT3_models_", run_id))
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(
  run_id = character(),
  model_id = character(),
  artifact_type = character(),
  file = character(),
  stringsAsFactors = FALSE
)

add_manifest <- function(model_id, artifact_type, file) {
  manifest <<- rbind(
    manifest,
    data.frame(
      run_id = run_id,
      model_id = model_id,
      artifact_type = artifact_type,
      file = file,
      stringsAsFactors = FALSE
    )
  )
}

safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  sub("^_|_$", "", x)
}

model_dir <- function(model_id) {
  d <- file.path(out_base, model_id)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

write_lines_safe <- function(lines, path) {
  writeLines(lines, con = path)
  message("Wrote: ", normalizePath(path))
}

write_csv_safe <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
  message("Wrote: ", normalizePath(path))
}

# ============================================================
# 0a) Preconditions / sanity
# ============================================================

stopifnot(exists("gdf_ct3"), exists("coords_ct3"), exists("coords_resid_ct3"))
stopifnot("specimen" %in% names(gdf_ct3))
stopifnot(identical(dimnames(coords_ct3)[[3]], gdf_ct3$specimen))
stopifnot(identical(dimnames(coords_resid_ct3), dimnames(coords_ct3)))

if (!("size_for_allometry" %in% names(gdf_ct3))) {
  stop("gdf_ct3 must contain size_for_allometry (set in R/methods/01_specimen_sampling_study_design/01_build_metadata.R).")
}

# Ensure CT3 grouping exists
if (!("group" %in% names(gdf_ct3))) {
  # Try to build it from year if needed
  if (!("year" %in% names(gdf_ct3))) stop("gdf_ct3 must contain `group` or `year`.")
  gdf_ct3 <- gdf_ct3 %>%
    mutate(group = case_when(
      year == 1950 ~ "CT_1950",
      year == 1956 ~ "CT_1956",
      year == 1970 ~ "CT_1970",
      TRUE ~ as.character(year)
    ))
}

# stable group ordering
gdf_ct3$group <- factor(gdf_ct3$group, levels = c("CT_1950", "CT_1956", "CT_1970"))

gdf_ct3$size_for_allometry <- as.numeric(gdf_ct3$size_for_allometry)

size_label_txt <- if ("size_label" %in% names(gdf_ct3)) unique(gdf_ct3$size_label)[1] else "size_for_allometry"

# ============================================================
# Helper: run procD.lm (geomorph version expects f1=)
# ============================================================

run_procD <- function(coords_arr, rhs, gdf_sub, iter = 999, RRPP = TRUE) {
  stopifnot(is.character(rhs), length(rhs) == 1)
  f <- stats::as.formula(paste("coords_arr ~", rhs))
  environment(f) <- environment()
  geomorph::procD.lm(f1 = f, data = gdf_sub, iter = iter, RRPP = RRPP)
}

# ============================================================
# Helper: write model card
# ============================================================

write_model_card <- function(model_id, coords_type, formula_text, gdf_sub,
                             size_covariate = NULL, size_label = NULL, notes = NULL) {
  d <- model_dir(model_id)
  path <- file.path(d, paste0(model_id, "_MODEL_CARD.txt"))
  
  lines <- c(
    paste0("run_id: ", run_id),
    paste0("model_id: ", model_id),
    paste0("coords_type: ", coords_type),
    paste0("formula: ", formula_text),
    "",
    paste0("n: ", nrow(gdf_sub)),
    "group_counts:",
    capture.output(print(table(gdf_sub$group, useNA = "ifany")))
  )
  
  if (!is.null(size_covariate)) lines <- c(lines, "", paste0("size_covariate: ", size_covariate))
  if (!is.null(size_label))     lines <- c(lines, paste0("size_label: ", size_label))
  if (!is.null(notes))          lines <- c(lines, "", "notes:", notes)
  
  lines <- c(lines, "", "sessionInfo():", capture.output(sessionInfo()))
  
  write_lines_safe(lines, path)
  add_manifest(model_id, "model_card_txt", path)
  invisible(path)
}

# ============================================================
# Helper: save fit + anova summary
# ============================================================

save_fit_bundle <- function(model_id, fit) {
  d <- model_dir(model_id)
  
  txt_path <- file.path(d, paste0(model_id, "_summary.txt"))
  lines <- c(
    paste0("run_id: ", run_id),
    paste0("model_id: ", model_id),
    "",
    "===== fit print =====",
    capture.output(print(fit)),
    "",
    "===== anova(fit) =====",
    capture.output(anova(fit))
  )
  write_lines_safe(lines, txt_path)
  add_manifest(model_id, "summary_txt", txt_path)
  
  rds_path <- file.path(d, paste0(model_id, "_fit.rds"))
  saveRDS(fit, rds_path)
  add_manifest(model_id, "fit_rds", rds_path)
  
  invisible(list(txt = txt_path, rds = rds_path))
}

# ============================================================
# Helper: mean shapes + simple Procrustes mean distance
# ============================================================

group_mean_shapes <- function(coords_arr, group_vec) {
  group_vec <- factor(group_vec)
  levs <- levels(group_vec)
  out <- setNames(vector("list", length(levs)), levs)
  
  for (g in levs) {
    idx <- which(group_vec == g)
    M <- apply(coords_arr[, , idx, drop = FALSE], c(1, 2), mean)
    rownames(M) <- dimnames(coords_arr)[[1]]
    colnames(M) <- c("x", "y")
    out[[g]] <- M
  }
  out
}

procdist_mean <- function(M1, M2) sqrt(sum((M1 - M2)^2))

write_mean_shape_distance_tables <- function(model_id, coords_arr, gdf_sub) {
  d <- model_dir(model_id)
  
  means <- group_mean_shapes(coords_arr, gdf_sub$group)
  grand <- geomorph::mshape(coords_arr)
  
  dist_to_grand <- data.frame(
    group = names(means),
    dist_to_grand = vapply(means, function(M) procdist_mean(M, grand), numeric(1)),
    stringsAsFactors = FALSE
  )
  csv1 <- file.path(d, paste0(model_id, "_meanShapeDist_toGrand.csv"))
  write_csv_safe(dist_to_grand, csv1)
  add_manifest(model_id, "csv_meanShapeDist_toGrand", csv1)
  
  levs <- names(means)
  pairs <- utils::combn(levs, 2)
  df_pairs <- data.frame(
    grp1 = pairs[1, ],
    grp2 = pairs[2, ],
    dist = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(df_pairs))) {
    df_pairs$dist[i] <- procdist_mean(means[[df_pairs$grp1[i]]], means[[df_pairs$grp2[i]]])
  }
  csv2 <- file.path(d, paste0(model_id, "_meanShapeDist_pairwise.csv"))
  write_csv_safe(df_pairs, csv2)
  add_manifest(model_id, "csv_meanShapeDist_pairwise", csv2)
  
  invisible(list(to_grand = csv1, pairwise = csv2))
}

# ============================================================
# Helper: morphol.disparity -> df
# ============================================================

disp_to_df <- function(disp_obj) {
  if (is.list(disp_obj) && !is.null(disp_obj$Procrustes.var)) {
    pv <- disp_obj$Procrustes.var
    return(data.frame(
      group = names(pv),
      Procrustes.var = as.numeric(pv),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(group = NA_character_, Procrustes.var = NA_real_, stringsAsFactors = FALSE)
}

# ============================================================
# Helper: morphol.disparity (ROBUST EVAL)
#   - assigns coords array into .GlobalEnv temporarily
# ============================================================

write_disparity_outputs <- function(model_id, coords_arr, gdf_sub, iter = 999) {
  d <- model_dir(model_id)
  
  tmp_name <- paste0("._tmp_coords_", safe_id(model_id), "_", as.integer(Sys.time()))
  assign(tmp_name, coords_arr, envir = .GlobalEnv)
  on.exit({
    if (exists(tmp_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = tmp_name, envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  f <- stats::as.formula(paste(tmp_name, "~ group"))
  
  disp <- geomorph::morphol.disparity(
    f1     = f,
    groups = gdf_sub$group,
    data   = gdf_sub,
    iter   = iter
  )
  
  txt_path <- file.path(d, paste0(model_id, "_morpholDisparity_print.txt"))
  write_lines_safe(capture.output(print(disp)), txt_path)
  add_manifest(model_id, "disparity_txt", txt_path)
  
  df <- disp_to_df(disp)
  csv_path <- file.path(d, paste0(model_id, "_morpholDisparity_byGroup.csv"))
  write_csv_safe(df, csv_path)
  add_manifest(model_id, "disparity_csv", csv_path)
  
  rds_path <- file.path(d, paste0(model_id, "_morpholDisparity.rds"))
  saveRDS(disp, rds_path)
  add_manifest(model_id, "disparity_rds", rds_path)
  
  invisible(list(txt = txt_path, csv = csv_path, rds = rds_path))
}

# ============================================================
# Helper: RRPP pairwise group tests (MEAN SHAPE differences)
# ============================================================

write_pairwise_rrpp <- function(model_id, fit, groups) {
  d <- model_dir(model_id)
  
  pw <- RRPP::pairwise(fit, groups = groups)
  
  rds_path <- file.path(d, paste0(model_id, "_pairwise_group_RRPP.rds"))
  saveRDS(pw, rds_path)
  add_manifest(model_id, "pairwise_rds", rds_path)
  
  txt_path <- file.path(d, paste0(model_id, "_pairwise_group_RRPP.txt"))
  sink(txt_path)
  on.exit(sink(), add = TRUE)
  cat("RRPP::pairwise output (print):\n\n")
  print(pw)
  cat("\n\nRRPP::pairwise output (summary):\n\n")
  s <- summary(pw, stat.table = TRUE, test.type = "dist")
  print(s)
  sink()
  add_manifest(model_id, "pairwise_txt", txt_path)
  
  tab <- NULL
  candidates <- list(
    s$summary.table, s$stat.table, s$tables, s$pairwise.tables, s$pairwise.table, s$pairwise,
    pw$summary.table, pw$stat.table, pw$tables
  )
  
  for (obj in candidates) {
    if (is.null(obj)) next
    if (is.data.frame(obj)) { tab <- obj; break }
    if (is.list(obj)) {
      for (x in obj) {
        if (is.data.frame(x)) { tab <- x; break }
      }
      if (!is.null(tab)) break
    }
  }
  
  if (!is.null(tab)) {
    csv_path <- file.path(d, paste0(model_id, "_pairwise_group_RRPP_table.csv"))
    utils::write.csv(tab, csv_path, row.names = FALSE)
    message("Wrote: ", normalizePath(csv_path))
    add_manifest(model_id, "pairwise_csv", csv_path)
  } else {
    debug_path <- file.path(d, paste0(model_id, "_pairwise_group_RRPP_STR.txt"))
    write_lines_safe(c(
      "Could not auto-extract a data.frame stat table from summary(pw).",
      "Here is str(pw):",
      capture.output(str(pw, max.level = 3)),
      "",
      "Here is str(summary(pw)):",
      capture.output(str(s, max.level = 3))
    ), debug_path)
    add_manifest(model_id, "pairwise_str_txt", debug_path)
  }
  
  invisible(list(pw = pw, summary_obj = s))
}

# ============================================================
# SIZE BUNDLE (CT3)
# ============================================================

model_id_size <- safe_id("SIZE_summaries_CT3")
write_model_card(
  model_id = model_id_size,
  coords_type = "RAW",
  formula_text = "N/A (descriptive size summaries / plots)",
  gdf_sub = gdf_ct3,
  size_covariate = "logCsize",
  size_label = "logCsize / centroid size from GPA",
  notes = c("Contains descriptive size summaries and a size distribution plot; not a procD model.")
)

d_size <- model_dir(model_id_size)

size_summary <- gdf_ct3 %>%
  group_by(group) %>%
  summarize(
    n = n(),
    mean_logCsize = mean(logCsize, na.rm = TRUE),
    sd_logCsize   = sd(logCsize, na.rm = TRUE),
    mean_Csize    = mean(Csize, na.rm = TRUE),
    sd_Csize      = sd(Csize, na.rm = TRUE),
    mean_SL_mm    = mean(SL_mm, na.rm = TRUE),
    sd_SL_mm      = sd(SL_mm, na.rm = TRUE),
    .groups = "drop"
  )

csv_size <- file.path(d_size, paste0(model_id_size, "_size_summary_by_group.csv"))
write_csv_safe(size_summary, csv_size)
add_manifest(model_id_size, "csv_size_summary", csv_size)

p_size <- ggplot(gdf_ct3, aes(group, logCsize, fill = group)) +
  geom_boxplot(alpha = 0.25, outlier.alpha = 0.6) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 1.6) +
  labs(
    title = "CT timeseries centroid size distribution by group",
    x = "Group",
    y = "log(Centroid size)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

png_size <- file.path(d_size, paste0(model_id_size, "_logCsize_by_group.png"))
ggsave(png_size, p_size, width = 7.5, height = 4.75, dpi = 300)
add_manifest(model_id_size, "figure_png_size", png_size)

# ============================================================
# MODELS (CT3)
# ============================================================

# MODEL 1: RAW shape ~ size
model_id_m1_raw <- safe_id("RAW_Model1_shape_by_size_CT3")
write_model_card(
  model_id = model_id_m1_raw,
  coords_type = "RAW",
  formula_text = "coords_ct3 ~ size_for_allometry",
  gdf_sub = gdf_ct3,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt
)
fit_raw_m1 <- run_procD(coords_ct3, "size_for_allometry", gdf_ct3, iter = 999)
save_fit_bundle(model_id_m1_raw, fit_raw_m1)

# MODEL 2: RAW shape ~ group
model_id_m2_raw <- safe_id("RAW_Model2_shape_by_group_CT3")
write_model_card(
  model_id = model_id_m2_raw,
  coords_type = "RAW",
  formula_text = "coords_ct3 ~ group",
  gdf_sub = gdf_ct3,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt,
  notes = c("Includes RRPP pairwise mean-shape comparisons, mean-shape distance tables, and morphol.disparity.")
)
fit_raw_m2 <- run_procD(coords_ct3, "group", gdf_ct3, iter = 999)
save_fit_bundle(model_id_m2_raw, fit_raw_m2)
write_pairwise_rrpp(model_id_m2_raw, fit_raw_m2, gdf_ct3$group)
write_mean_shape_distance_tables(model_id_m2_raw, coords_ct3, gdf_ct3)
write_disparity_outputs(model_id_m2_raw, coords_ct3, gdf_ct3)

# MODEL 3: RAW shape ~ size * group
model_id_m3_raw <- safe_id("RAW_Model3_shape_by_sizeXgroup_CT3")
write_model_card(
  model_id = model_id_m3_raw,
  coords_type = "RAW",
  formula_text = "coords_ct3 ~ size_for_allometry * group",
  gdf_sub = gdf_ct3,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt
)
fit_raw_m3 <- run_procD(coords_ct3, "size_for_allometry * group", gdf_ct3, iter = 999)
save_fit_bundle(model_id_m3_raw, fit_raw_m3)

# MODEL 2 (RESID): residual shape ~ group
model_id_m2_res <- safe_id("RESID_Model2_shape_by_group_CT3")
write_model_card(
  model_id = model_id_m2_res,
  coords_type = "RESID",
  formula_text = "coords_resid_ct3 ~ group",
  gdf_sub = gdf_ct3,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt,
  notes = c("Includes RRPP pairwise mean-shape comparisons, mean-shape distance tables, and morphol.disparity (on residual shapes).")
)
fit_res_m2 <- run_procD(coords_resid_ct3, "group", gdf_ct3, iter = 999)
save_fit_bundle(model_id_m2_res, fit_res_m2)
write_pairwise_rrpp(model_id_m2_res, fit_res_m2, gdf_ct3$group)
write_mean_shape_distance_tables(model_id_m2_res, coords_resid_ct3, gdf_ct3)
write_disparity_outputs(model_id_m2_res, coords_resid_ct3, gdf_ct3)


# ============================================================
# SIZE DIAGNOSTIC (CT3):
# Restrict to overlapping size range and re-run ANCOVA
#   - Filter: 2.9 < logCsize < 3.35
#   - Refit: shape ~ size_for_allometry * group
# Saves to: Outputs/<run_id>/SIZE_diagnostic/
# ============================================================

# ---- settings ----
DIAG_LO <- 2.9
DIAG_HI <- 3.35
DIAG_ITER <- 999

# ---- helper: subset coords by gdf row order (specimens) ----
subset_coords_by_specimens <- function(coords_arr, gdf_sub) {
  stopifnot("specimen" %in% names(gdf_sub))
  keep <- as.character(gdf_sub$specimen)
  stopifnot(all(keep %in% dimnames(coords_arr)[[3]]))
  coords_arr[, , keep, drop = FALSE]
}

# ---- create diagnostic folder (inside this run) ----
diag_model_id <- safe_id("SIZE_diagnostic_overlapRange_CT3")
diag_dir <- model_dir(diag_model_id)

# ---- filter gdf to overlapping size window ----
gdf_ct3_diag <- gdf_ct3 %>%
  dplyr::filter(
    !is.na(logCsize),
    logCsize > DIAG_LO,
    logCsize < DIAG_HI
  ) %>%
  droplevels()

# ---- subset coords + resid coords to same specimens ----
coords_ct3_diag       <- subset_coords_by_specimens(coords_ct3, gdf_ct3_diag)
coords_resid_ct3_diag <- subset_coords_by_specimens(coords_resid_ct3, gdf_ct3_diag)

# ---- sanity ----
stopifnot(identical(dimnames(coords_ct3_diag)[[3]], gdf_ct3_diag$specimen))
stopifnot(identical(dimnames(coords_resid_ct3_diag), dimnames(coords_ct3_diag)))

# ---- write a short diagnostic card (counts + size range by group) ----
diag_card <- file.path(diag_dir, paste0(diag_model_id, "_DIAGNOSTIC_CARD.txt"))
diag_lines <- c(
  paste0("run_id: ", run_id),
  paste0("model_id: ", diag_model_id),
  "",
  "Purpose: test whether size*group interaction is driven by non-overlapping size distributions.",
  paste0("Filter applied: logCsize > ", DIAG_LO, " and logCsize < ", DIAG_HI),
  "",
  paste0("n (before filter): ", nrow(gdf_ct3)),
  paste0("n (after filter):  ", nrow(gdf_ct3_diag)),
  "",
  "group_counts (after filter):",
  capture.output(print(table(gdf_ct3_diag$group, useNA = "ifany"))),
  "",
  "logCsize summary by group (after filter):",
  capture.output(
    print(
      gdf_ct3_diag %>%
        group_by(group) %>%
        summarize(
          n = n(),
          min_logCsize = min(logCsize, na.rm = TRUE),
          max_logCsize = max(logCsize, na.rm = TRUE),
          mean_logCsize = mean(logCsize, na.rm = TRUE),
          sd_logCsize = sd(logCsize, na.rm = TRUE),
          .groups = "drop"
        )
    )
  ),
  "",
  "NOTE: Refit uses the same RHS as RAW_Model3_shape_by_sizeXgroup_CT3, but on filtered specimens."
)
write_lines_safe(diag_lines, diag_card)
add_manifest(diag_model_id, "diagnostic_card_txt", diag_card)

# ---- optional: save the filtered gdf as CSV ----
diag_csv <- file.path(diag_dir, paste0(diag_model_id, "_gdf_filtered.csv"))
write_csv_safe(gdf_ct3_diag, diag_csv)
add_manifest(diag_model_id, "gdf_filtered_csv", diag_csv)

# ---- optional: size distribution plot for filtered set ----
p_diag <- ggplot(gdf_ct3_diag, aes(group, logCsize, fill = group)) +
  geom_boxplot(alpha = 0.25, outlier.alpha = 0.6) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 1.6) +
  labs(
    title = paste0("CT3 size diagnostic (filtered): ", DIAG_LO, " < logCsize < ", DIAG_HI),
    x = "Group",
    y = "log(Centroid size)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

png_diag <- file.path(diag_dir, paste0(diag_model_id, "_logCsize_by_group_filtered.png"))
ggsave(png_diag, p_diag, width = 7.5, height = 4.75, dpi = 300)
add_manifest(diag_model_id, "figure_png_size_filtered", png_diag)

# ---- refit ANCOVA on filtered RAW shapes ----
diag_raw_m3_id <- safe_id("RAW_Model3_shape_by_sizeXgroup_CT3_overlapRange")
write_model_card(
  model_id = diag_raw_m3_id,
  coords_type = "RAW",
  formula_text = "coords_ct3_diag ~ size_for_allometry * group  (FILTERED: 2.9<logCsize<3.35)",
  gdf_sub = gdf_ct3_diag,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt,
  notes = c(
    "Size diagnostic (gold-standard check): restrict to overlapping size range across groups.",
    paste0("Filter: ", DIAG_LO, " < logCsize < ", DIAG_HI),
    "Interpretation: if size:group interaction weakens/disappears, interaction in full dataset may be driven by non-overlapping size distributions."
  )
)
fit_raw_m3_diag <- run_procD(coords_ct3_diag, "size_for_allometry * group", gdf_ct3_diag, iter = DIAG_ITER)
save_fit_bundle(diag_raw_m3_id, fit_raw_m3_diag)

# ---- optional: refit on residual shapes too (same filter) ----
diag_res_m3_id <- safe_id("RESID_Model3_shape_by_sizeXgroup_CT3_overlapRange")
write_model_card(
  model_id = diag_res_m3_id,
  coords_type = "RESID",
  formula_text = "coords_resid_ct3_diag ~ size_for_allometry * group  (FILTERED: 2.9<logCsize<3.35)",
  gdf_sub = gdf_ct3_diag,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt,
  notes = c(
    "Same diagnostic as RAW, but applied to size-corrected residual shapes (if you want to confirm the pattern there too).",
    paste0("Filter: ", DIAG_LO, " < logCsize < ", DIAG_HI)
  )
)
fit_res_m3_diag <- run_procD(coords_resid_ct3_diag, "size_for_allometry * group", gdf_ct3_diag, iter = DIAG_ITER)
save_fit_bundle(diag_res_m3_id, fit_res_m3_diag)

message("SIZE DIAGNOSTIC COMPLETE. Folder: ", normalizePath(diag_dir))




# ============================================================
# MANIFEST
# ============================================================

manifest_path <- file.path(out_base, "MANIFEST.csv")
write_csv_safe(manifest, manifest_path)
message("DONE. Outputs in: ", normalizePath(out_base))